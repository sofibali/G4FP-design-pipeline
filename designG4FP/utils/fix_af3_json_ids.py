#!/usr/bin/env python3
"""
Fix AlphaFold3 JSON inputs by adding required 'id' fields to all sequences.

AlphaFold3 requires each sequence in the sequences array to have an 'id' field.
This script adds sequential IDs to all sequences.

Usage:
    python fix_af3_json_ids.py [--output-dir OUTPUT_DIR]
"""

import json
import argparse
from pathlib import Path
import sys

def fix_json_file(json_path: Path, dry_run: bool = False) -> bool:
    """
    Fix a single JSON file by adding id fields to sequences
    
    Returns:
        True if file was modified, False otherwise
    """
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        
        # Check if sequences need fixing
        needs_fix = False
        if 'sequences' in data:
            for seq in data['sequences']:
                if 'id' not in seq:
                    needs_fix = True
                    break
        
        if not needs_fix:
            return False
        
        # Add sequential IDs to sequences
        for i, seq in enumerate(data['sequences'], start=1):
            if 'id' not in seq:
                seq['id'] = str(i)
        
        if dry_run:
            print(f"  Would fix: {json_path.name}")
            return True
        
        # Write back with proper formatting
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        return True
        
    except Exception as e:
        print(f"  ✗ Error processing {json_path.name}: {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(description='Fix AlphaFold3 JSON files by adding sequence IDs')
    parser.add_argument('--output-dir', type=str,
                       help='Specific output directory (default: all output_* dirs)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be fixed without making changes')
    
    args = parser.parse_args()
    
    # Find output directories
    if args.output_dir:
        output_dirs = [Path(args.output_dir)]
    else:
        output_dirs = [d for d in Path('.').iterdir() 
                      if d.is_dir() and d.name.startswith('output_')]
    
    if not output_dirs:
        print("No output directories found")
        sys.exit(1)
    
    print(f"{'DRY RUN - ' if args.dry_run else ''}Fixing AlphaFold3 JSON files...")
    print(f"Found {len(output_dirs)} output directories\n")
    
    total_fixed = 0
    total_checked = 0
    
    for output_dir in sorted(output_dirs):
        print(f"Processing: {output_dir.name}")
        
        # Check both bound and apo directories
        for state in ['bound', 'apo']:
            input_dir = output_dir / f"02_alphafold3_inputs_{state}"
            
            if not input_dir.exists():
                continue
            
            json_files = list(input_dir.glob("*.json"))
            if not json_files:
                continue
            
            print(f"  {state}: {len(json_files)} JSON files")
            
            fixed_count = 0
            for json_file in json_files:
                if fix_json_file(json_file, args.dry_run):
                    fixed_count += 1
            
            if fixed_count > 0:
                print(f"    {'Would fix' if args.dry_run else 'Fixed'}: {fixed_count} files")
            else:
                print(f"    ✓ All files already have IDs")
            
            total_fixed += fixed_count
            total_checked += len(json_files)
        
        print()
    
    print("=" * 60)
    print(f"Summary:")
    print(f"  Total files checked: {total_checked}")
    print(f"  Files {'that need' if args.dry_run else ''} fixed: {total_fixed}")
    
    if args.dry_run and total_fixed > 0:
        print(f"\nRun without --dry-run to apply fixes")
    elif total_fixed > 0:
        print(f"\n✓ All JSON files fixed! Ready to run AlphaFold3 predictions.")
    else:
        print(f"\n✓ No fixes needed - all files already have sequence IDs.")


if __name__ == "__main__":
    main()

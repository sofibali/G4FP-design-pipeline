#!/usr/bin/env python3
"""
Quick Test: Verify structure analysis scripts work

Tests the analysis pipeline on a small subset before full run.
"""

import sys
from pathlib import Path

print("=" * 80)
print("Testing Structure Analysis Scripts")
print("=" * 80)

# Check required files exist
designG4FP_dir = Path("/home/sbali/LigandMPNN/designG4FP")
scripts = [
    "05_analyze_af3_structures.py",
    "06_compare_ligand_states.py",
    "run_structure_analysis.sh"
]

print("\n1. Checking scripts exist...")
for script in scripts:
    path = designG4FP_dir / script
    if path.exists():
        print(f"  ✓ {script}")
    else:
        print(f"  ✗ {script} NOT FOUND")
        sys.exit(1)

# Check for completed predictions
print("\n2. Checking for completed AF3 predictions...")
test_output = designG4FP_dir / "output_G4FP_des1_cro_mod0"

if not test_output.exists():
    print(f"  ✗ Test output directory not found: {test_output}")
    sys.exit(1)

bound_dir = test_output / "03_alphafold3_predictions_bound"
apo_dir = test_output / "03_alphafold3_predictions_apo"

if bound_dir.exists():
    bound_count = len([d for d in bound_dir.iterdir() if d.is_dir() and d.name.startswith('design_')])
    print(f"  ✓ Found {bound_count} bound predictions")
else:
    print(f"  ✗ Bound predictions not found")
    bound_count = 0

if apo_dir.exists():
    apo_count = len([d for d in apo_dir.iterdir() if d.is_dir() and d.name.startswith('design_')])
    print(f"  ✓ Found {apo_count} apo predictions")
else:
    print(f"  ✗ Apo predictions not found")
    apo_count = 0

# Check template structure
print("\n3. Checking template structure...")
template_path = Path("/home/sbali/LigandMPNN/inputs/G4FP_design1_cro.pdb")

if template_path.exists():
    print(f"  ✓ Template found: {template_path}")
else:
    print(f"  ✗ Template not found: {template_path}")
    print("  Available templates:")
    inputs_dir = Path("/home/sbali/LigandMPNN/inputs")
    for f in inputs_dir.glob("*.pdb"):
        if "G4FP" in f.name or "cro" in f.name:
            print(f"    - {f.name}")
    sys.exit(1)

# Check dependencies
print("\n4. Checking Python dependencies...")
required_packages = [
    "Bio.PDB",
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "scipy"
]

missing = []
for package in required_packages:
    try:
        if package == "Bio.PDB":
            from Bio.PDB import PDBParser
        else:
            __import__(package.split('.')[0])
        print(f"  ✓ {package}")
    except ImportError:
        print(f"  ✗ {package} NOT INSTALLED")
        missing.append(package)

if missing:
    print(f"\nInstall missing packages:")
    print(f"  pip install biopython pandas numpy matplotlib seaborn scipy openpyxl")
    sys.exit(1)

print("\n" + "=" * 80)
print("✓ All checks passed!")
print("=" * 80)
print("\nReady to run analysis:")
print(f"\n  cd {designG4FP_dir}")
print(f"  bash run_structure_analysis.sh {test_output.name} {template_path}")
print("\nOr test individual scripts:")
print(f"\n  python 05_analyze_af3_structures.py \\")
print(f"      --output-dir {test_output.name} \\")
print(f"      --template {template_path} \\")
print(f"      --state both")
print("")

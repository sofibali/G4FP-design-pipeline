#!/usr/bin/env python3
"""
Build a PyMOL session (.pse) with the top-ranked G4FP designs.

Loads:
  - Template PDB (input structure, white)
  - Template AF3 holo + apo best models (if available)
  - Top N designs: holo (warm colors) and apo (cool colors)
  - Colors by pLDDT (B-factor spectrum)

Usage:
    source /programs/sbgrid.shrc
    pymol -cq build_pymol_session.py
    # or: pymol -cq build_pymol_session.py -- --top 10 --output top10.pse
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict

# Parse args after "--"
top_n = 10
output_name = "07_results/07_top_designs.pse"
args = sys.argv[1:]
if "--" in sys.argv:
    args = sys.argv[sys.argv.index("--") + 1:]
i = 0
while i < len(args):
    if args[i] == "--top" and i + 1 < len(args):
        top_n = int(args[i + 1])
        i += 2
    elif args[i] == "--output" and i + 1 < len(args):
        output_name = args[i + 1]
        i += 2
    else:
        i += 1

# ── Setup ──
try:
    from pymol import cmd, stored
except ImportError:
    print("ERROR: Must run inside PyMOL (pymol -cq build_pymol_session.py)")
    sys.exit(1)

import csv

# __file__ may not resolve correctly when run via pymol -cq; use cwd as fallback
try:
    _script_path = Path(__file__).resolve()
    base_dir = _script_path.parent
    if not (base_dir / "07_results").exists():
        base_dir = Path.cwd()
except NameError:
    base_dir = Path.cwd()

# ── Load top candidates ──
candidates_csv = base_dir / "07_results" / "07_final_candidates.csv"
if not candidates_csv.exists():
    print(f"ERROR: {candidates_csv} not found. Run 07_aggregate_and_rank.py first.")
    sys.exit(1)

with open(candidates_csv) as f:
    reader = csv.DictReader(f)
    candidates = list(reader)

top = candidates[:top_n]
print(f"\nBuilding PyMOL session with top {len(top)} designs...")

# ── Discover structure files ──

def find_design_dir(predictions_dir, seq_id, state):
    """Find the AF3 result directory for a design (handles timestamped dirs)."""
    pattern = re.compile(
        rf'^design_{seq_id}_{re.escape(state)}(?:_(\d{{8}}_\d{{6}}))?$'
    )
    candidates = []
    if not predictions_dir.exists():
        return None
    for d in predictions_dir.iterdir():
        if not d.is_dir():
            continue
        m = pattern.match(d.name)
        if m:
            has_model = bool(list(d.glob("*_model.cif")))
            if has_model:
                ts = m.group(1)
                candidates.append((ts is not None, ts or "", d))
    if not candidates:
        return None
    # Prefer non-timestamped, then latest timestamp
    non_ts = [e for e in candidates if not e[0]]
    if non_ts:
        return non_ts[0][2]
    candidates.sort(key=lambda e: e[1], reverse=True)
    return candidates[0][2]


def get_best_model(design_dir):
    """Get the best (top-level) model CIF from a design directory."""
    cifs = list(design_dir.glob("*_model.cif"))
    if cifs:
        return cifs[0]
    return None


# ── Color scheme ──
# Warm palette for holo (oranges/reds), cool for apo (blues/teals)
holo_colors = [
    [1.00, 0.60, 0.20],  # orange
    [0.95, 0.45, 0.15],
    [0.90, 0.35, 0.10],
    [0.85, 0.28, 0.08],
    [0.80, 0.20, 0.05],
    [0.95, 0.55, 0.25],
    [0.88, 0.40, 0.12],
    [0.82, 0.32, 0.09],
    [0.78, 0.25, 0.06],
    [0.75, 0.18, 0.04],
    [1.00, 0.50, 0.18],
    [0.92, 0.42, 0.14],
    [0.86, 0.34, 0.10],
    [0.80, 0.26, 0.07],
    [0.76, 0.20, 0.05],
    [0.98, 0.58, 0.22],
    [0.90, 0.38, 0.13],
    [0.84, 0.30, 0.09],
    [0.79, 0.23, 0.06],
    [0.74, 0.17, 0.04],
]

apo_colors = [
    [0.30, 0.60, 0.90],  # blue
    [0.25, 0.50, 0.85],
    [0.20, 0.45, 0.80],
    [0.15, 0.40, 0.75],
    [0.10, 0.35, 0.70],
    [0.35, 0.55, 0.88],
    [0.28, 0.48, 0.82],
    [0.22, 0.42, 0.78],
    [0.17, 0.38, 0.73],
    [0.12, 0.33, 0.68],
    [0.32, 0.58, 0.92],
    [0.26, 0.52, 0.86],
    [0.21, 0.46, 0.80],
    [0.16, 0.40, 0.76],
    [0.11, 0.34, 0.71],
    [0.33, 0.56, 0.89],
    [0.27, 0.49, 0.83],
    [0.22, 0.43, 0.79],
    [0.17, 0.37, 0.74],
    [0.12, 0.32, 0.69],
]

# ── Reinitialize PyMOL ──
cmd.reinitialize()
cmd.set("cif_use_auth", 0)

# ── 1. Load template input PDB ──
template_name = top[0].get("template", "G4FP_des1_cro_mod0")
output_dir_name = f"output_{template_name}" if not template_name.startswith("output_") else template_name
template_name_clean = template_name.replace("output_", "")

template_pdb = base_dir / "inputs" / f"{template_name_clean}.pdb"
if template_pdb.exists():
    obj_name = "template_input"
    cmd.load(str(template_pdb), obj_name)
    cmd.color("white", obj_name)
    cmd.show("cartoon", obj_name)
    cmd.set("cartoon_transparency", 0.5, obj_name)
    print(f"  Loaded template PDB: {template_pdb.name}")
else:
    print(f"  WARNING: Template PDB not found: {template_pdb}")

# ── 2. Load template AF3 structures ──
template_af3_holo = base_dir / "template_af3_outputs" / "bound"  # disk dir is "bound"
template_af3_apo = base_dir / "template_af3_outputs" / "apo"

for display_state, disk_state, af3_dir, color in [
    ("holo", "bound", template_af3_holo, [0.4, 0.8, 0.4]),   # green
    ("apo", "apo", template_af3_apo, [0.7, 0.4, 0.7]),       # purple
]:
    # Find matching template dir (lowercase variant)
    tpl_pattern = f"template_{template_name_clean.lower()}_{disk_state}"
    found = None
    if af3_dir.exists():
        for d in af3_dir.iterdir():
            if d.is_dir() and d.name.lower().replace("_", "") == tpl_pattern.replace("_", ""):
                found = d
                break
            # Also try direct match
            if d.is_dir() and d.name.startswith(f"template_") and disk_state in d.name:
                inner_name = template_name_clean.lower().replace("_", "")
                dir_name = d.name.lower().replace("_", "").replace("template", "").replace(disk_state, "")
                if inner_name.replace("_", "") == dir_name:
                    found = d
                    break
        # Simpler: just glob for it
        if found is None:
            for d in af3_dir.iterdir():
                if d.is_dir() and template_name_clean.lower().replace("_", "") in d.name.lower().replace("_", ""):
                    found = d
                    break

    if found is not None:
        # Check for nested dir
        model_cif = list(found.glob("*_model.cif"))
        if not model_cif:
            for sub in found.iterdir():
                if sub.is_dir():
                    model_cif = list(sub.glob("*_model.cif"))
                    if model_cif:
                        break
        if model_cif:
            obj_name = f"template_af3_{display_state}"
            cmd.load(str(model_cif[0]), obj_name)
            cmd.color("green" if display_state == "holo" else "purpleblue", obj_name)
            cmd.show("cartoon", obj_name)
            cmd.set("cartoon_transparency", 0.3, obj_name)
            print(f"  Loaded template AF3 {display_state}: {model_cif[0].name}")

# ── 3. Load top N designs (holo + apo) ──
predictions_holo = base_dir / output_dir_name / "03_alphafold3_predictions_bound"
predictions_apo = base_dir / output_dir_name / "03_alphafold3_predictions_apo"

loaded_holo = []
loaded_apo = []

for idx, row in enumerate(top):
    seq_id = int(row["seq_id"])
    rank = int(row["final_rank"])
    fitness = float(row["fitness_score"])

    for display_state, disk_state, pred_dir, colors, loaded_list in [
        ("holo", "bound", predictions_holo, holo_colors, loaded_holo),
        ("apo", "apo", predictions_apo, apo_colors, loaded_apo),
    ]:
        design_dir = find_design_dir(pred_dir, seq_id, disk_state)
        if design_dir is None:
            print(f"  WARNING: design_{seq_id}_{display_state} not found")
            continue

        model_cif = get_best_model(design_dir)
        if model_cif is None:
            continue

        obj_name = f"rank{rank:02d}_d{seq_id:04d}_{display_state}"
        cmd.load(str(model_cif), obj_name)

        # Color by state
        c = colors[idx % len(colors)]
        cmd.set_color(f"c_{obj_name}", c)
        cmd.color(f"c_{obj_name}", obj_name)
        cmd.show("cartoon", obj_name)

        loaded_list.append(obj_name)

    print(f"  Rank {rank}: design_{seq_id} "
          f"(fitness={fitness:.3f}, dpLDDT={float(row.get('mean_plddt_diff', 0)):.1f})")

# ── 4. Create groups ──
if loaded_holo:
    cmd.group("designs_holo", " ".join(loaded_holo))
if loaded_apo:
    cmd.group("designs_apo", " ".join(loaded_apo))

# Group templates
tpl_objs = []
if "template_input" in cmd.get_names():
    tpl_objs.append("template_input")
if "template_af3_holo" in cmd.get_names():
    tpl_objs.append("template_af3_holo")
if "template_af3_apo" in cmd.get_names():
    tpl_objs.append("template_af3_apo")
if tpl_objs:
    cmd.group("templates", " ".join(tpl_objs))

# ── 5. Align everything to template ──
ref_obj = "template_input" if "template_input" in cmd.get_names() else (
    loaded_holo[0] if loaded_holo else None)

if ref_obj:
    all_objs = cmd.get_names("objects")
    for obj in all_objs:
        if obj != ref_obj and obj not in ["templates", "designs_holo", "designs_apo"]:
            try:
                cmd.align(obj, ref_obj)
            except Exception:
                pass

# ── 6. Style settings ──
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("cartoon_oval_length", 1.2)
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.bg_color("white")

# Show chain A cartoon only, hide other chains for clarity
cmd.hide("everything")
cmd.show("cartoon", "chain A")

# Highlight chromophore region (residues 197-199) as sticks
cmd.select("chromophore", "resi 197-199 and chain A")
cmd.show("sticks", "chromophore")
cmd.color("yellow", "chromophore")
cmd.deselect()

# Start with apo hidden (user can toggle)
if loaded_apo:
    cmd.disable("designs_apo")

# Zoom to fit
cmd.zoom("chain A")
cmd.orient()

# ── 7. Save session ──
output_path = str(base_dir / output_name)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
cmd.save(output_path)

print(f"\nPyMOL session saved: {output_path}")
print(f"  {len(loaded_holo)} holo + {len(loaded_apo)} apo structures")
print(f"  Groups: templates, designs_holo, designs_apo")
print(f"  Apo group is disabled by default (toggle in PyMOL)")
print(f"\nOpen with: pymol {output_name}")

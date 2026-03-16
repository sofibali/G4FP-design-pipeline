#!/usr/bin/env python3
"""
Build a PyMOL session (.pse) with G4FP design candidates.

Loads:
  - Input template PDBs + 5 AF3 replicate structures per template
  - Top 3 designs per template (holo, 5 replicates each)
  - Top 10 global designs (holo, 5 replicates each), grouped by template
  - DNA (G-quadruplex) shown as cartoon
  - K+ ions shown as nb_spheres, colored by atom type
  - Chromophore (TYG 197-199) as a single selection group
  - G4 DNA as a single selection group

Usage:
    source /programs/sbgrid.shrc
    pymol -cq build_pymol_session.py
    pymol -cq build_pymol_session.py -- --output custom.pse
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict

# Parse args after "--"
output_name = "07_results/07_top_designs.pse"
N_REPLICATES = 5  # seed replicates to load per structure
TOP_PER_TEMPLATE = 3
TOP_GLOBAL = 10

args = sys.argv[1:]
if "--" in sys.argv:
    args = sys.argv[sys.argv.index("--") + 1:]
i = 0
while i < len(args):
    if args[i] == "--output" and i + 1 < len(args):
        output_name = args[i + 1]
        i += 2
    elif args[i] == "--replicates" and i + 1 < len(args):
        N_REPLICATES = int(args[i + 1])
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

try:
    _script_path = Path(__file__).resolve()
    base_dir = _script_path.parent
    if not (base_dir / "07_results").exists():
        base_dir = Path.cwd()
except NameError:
    base_dir = Path.cwd()

# ── tab10 colors matching matplotlib plots ──
TEMPLATE_COLORS = {
    "G4FP_des1_cro_mod0":   [0.1216, 0.4667, 0.7059],  # tab:blue
    "G4FP_des1_cro_mod4":   [1.0000, 0.4980, 0.0549],  # tab:orange
    "G4FP_r20_l1_cro_mod0": [0.1725, 0.6275, 0.1725],  # tab:green
    "G4FP_r20_l2_cro_mod2": [0.8392, 0.1529, 0.1569],  # tab:red
    "G4FP_r20_l2_cro_mod3": [0.5804, 0.4039, 0.7412],  # tab:purple
    "G4FP_r28_l1_cro_mod2": [0.5490, 0.3373, 0.2941],  # tab:brown
    "G4FP_r28_l2_cro_mod2": [0.8902, 0.4667, 0.7608],  # tab:pink
    "G4FP_r28_l2_cro_mod3": [0.4980, 0.4980, 0.4980],  # tab:gray
    "G4FP_r28_l3_cro_mod0": [0.7373, 0.7412, 0.1333],  # tab:olive
    "G4FP_r28_l3_cro_mod2": [0.0902, 0.7451, 0.8118],  # tab:cyan
}

# Chromophore is at CRO residue 197 in all templates → TYG at 197-199
CHROMOPHORE_RESI = "resi 197-199"


def discover_seed_dirs(design_dir):
    """Find seed-N_sample-M subdirectories, return sorted list of Paths."""
    pattern = re.compile(r'^seed-(\d+)_sample-(\d+)$')
    dirs = []
    for d in design_dir.iterdir():
        if d.is_dir() and pattern.match(d.name):
            dirs.append(d)
    if not dirs:
        # Check one level deeper (padded layout)
        for subdir in design_dir.iterdir():
            if subdir.is_dir():
                for d in subdir.iterdir():
                    if d.is_dir() and pattern.match(d.name):
                        dirs.append(d)
                if dirs:
                    break
    return sorted(dirs)[:N_REPLICATES]


def find_design_dir(predictions_dir, seq_id, state):
    """Find the AF3 output directory for a design."""
    if not predictions_dir.exists():
        return None
    # Try non-padded: design_N_state/
    non_padded = predictions_dir / f"design_{seq_id}_{state}"
    if non_padded.is_dir():
        return non_padded
    # Try padded: design_NNNN_state/
    padded = predictions_dir / f"design_{int(seq_id):04d}_{state}"
    if padded.is_dir():
        # Check for inner dir
        for sub in padded.iterdir():
            if sub.is_dir() and (list(sub.glob("seed-*")) or
                                  list(sub.glob("*_model.cif"))):
                return sub
        return padded
    return None


def find_template_af3_dir(tpl_name, state):
    """Find AF3 output dir for a template's bound/apo prediction."""
    disk_state = "bound" if state == "holo" else "apo"
    tpl_af3 = base_dir / "template_af3_outputs" / disk_state
    if not tpl_af3.exists():
        return None
    # Try various naming conventions
    for prefix in [f"template_{tpl_name}_{disk_state}",
                   f"template_{tpl_name.lower()}_{disk_state}"]:
        candidate = tpl_af3 / prefix
        if candidate.exists():
            # Check for inner dir with seeds
            if list(candidate.glob("seed-*")):
                return candidate
            for sub in candidate.iterdir():
                if sub.is_dir() and list(sub.glob("seed-*")):
                    return sub
    return None


# ── Load candidates ──
candidates_csv = base_dir / "07_results" / "07_final_candidates.csv"
if not candidates_csv.exists():
    print(f"ERROR: {candidates_csv} not found. Run 07_aggregate_and_rank.py first.")
    sys.exit(1)

with open(candidates_csv) as f:
    reader = csv.DictReader(f)
    all_candidates = list(reader)

# ── Select: top 3 per template + top 10 global ──
# Top 10 global (already sorted by fitness_score / final_rank)
top_global = all_candidates[:TOP_GLOBAL]

# Top 3 per template
by_template = defaultdict(list)
for row in all_candidates:
    by_template[row["template"]].append(row)
top_per_tpl = {}
for tpl, rows in by_template.items():
    top_per_tpl[tpl] = rows[:TOP_PER_TEMPLATE]

# Union of both sets (avoid duplicates)
selected = {}
for row in top_global:
    key = (row["template"], row["seq_id"])
    selected[key] = row
for tpl, rows in top_per_tpl.items():
    for row in rows:
        key = (row["template"], row["seq_id"])
        if key not in selected:
            selected[key] = row

templates_seen = sorted(set(row["template"] for row in selected.values()))
print(f"\nBuilding PyMOL session:")
print(f"  Top {TOP_GLOBAL} global + top {TOP_PER_TEMPLATE}/template "
      f"= {len(selected)} unique designs")
print(f"  {N_REPLICATES} AF3 replicates each")
print(f"  {len(templates_seen)} templates")

# ── Reinitialize PyMOL ──
cmd.reinitialize()
cmd.set("cif_use_auth", 1)

# Register template colors
for tpl, rgb in TEMPLATE_COLORS.items():
    cmd.set_color(f"tpl_{tpl}", rgb)

# ── 1. Load template structures ──
print("\nLoading templates...")
template_objs = []
tpl_replicate_objs = []

for tpl in templates_seen:
    rgb = TEMPLATE_COLORS.get(tpl, [0.7, 0.7, 0.7])
    color_name = f"tpl_{tpl}"
    short = tpl.replace("G4FP_", "")

    # 1a. Original input PDB
    template_pdb = base_dir / "inputs" / f"{tpl}.pdb"
    if template_pdb.exists():
        obj_name = f"tpl_input_{short}"
        cmd.load(str(template_pdb), obj_name)
        cmd.color(color_name, obj_name)
        cmd.show("cartoon", obj_name)
        cmd.set("cartoon_transparency", 0.5, obj_name)
        template_objs.append(obj_name)
        print(f"  Input PDB: {tpl}")

    # 1b. AF3 replicate structures (holo only, 5 replicates)
    af3_dir = find_template_af3_dir(tpl, "holo")
    if af3_dir is not None:
        seed_dirs = discover_seed_dirs(af3_dir)
        for j, sd in enumerate(seed_dirs):
            cif = sd / "model.cif"
            if not cif.exists():
                continue
            obj_name = f"tpl_af3_{short}_s{j}"
            cmd.load(str(cif), obj_name)
            cmd.color(color_name, obj_name)
            cmd.show("cartoon", f"{obj_name} and chain A")
            cmd.set("cartoon_transparency", 0.4, obj_name)
            tpl_replicate_objs.append(obj_name)
        if seed_dirs:
            print(f"  AF3 replicates: {tpl} ({len(seed_dirs)} seeds)")

# Group templates
all_tpl_objs = template_objs + tpl_replicate_objs
if all_tpl_objs:
    cmd.group("templates", " ".join(all_tpl_objs))
if template_objs:
    cmd.group("tpl_input_pdbs", " ".join(template_objs))
if tpl_replicate_objs:
    cmd.group("tpl_af3_replicates", " ".join(tpl_replicate_objs))

# ── 2. Load selected designs (holo, 5 replicates each) ──
print("\nLoading designs...")
all_design_objs = []
top10_objs = []
loaded_by_template = defaultdict(list)

# Track which are in top 10 global
top10_keys = set((row["template"], row["seq_id"]) for row in top_global)

for (tpl, seq_id), row in sorted(selected.items()):
    rank = int(row["final_rank"])
    fitness = float(row["fitness_score"])
    rgb = TEMPLATE_COLORS.get(tpl, [0.7, 0.7, 0.7])
    color_name = f"tpl_{tpl}"
    short_tpl = tpl.replace("G4FP_", "")
    output_dir = base_dir / f"output_{tpl}"

    # Only holo (bound) structures
    pred_dir = output_dir / "03_alphafold3_predictions_bound"
    design_dir = find_design_dir(pred_dir, seq_id, "bound")
    if design_dir is None:
        print(f"  SKIP: {tpl} design_{seq_id} (no predictions found)")
        continue

    # Load replicates
    seed_dirs = discover_seed_dirs(design_dir)
    if not seed_dirs:
        # Fall back to top-level model.cif
        top_cif = list(design_dir.glob("*_model.cif"))
        if top_cif:
            obj_name = f"r{rank:03d}_{short_tpl}_d{int(seq_id):04d}"
            cmd.load(str(top_cif[0]), obj_name)
            cmd.color(color_name, obj_name)
            all_design_objs.append(obj_name)
            loaded_by_template[tpl].append(obj_name)
            if (tpl, seq_id) in top10_keys:
                top10_objs.append(obj_name)
        continue

    for j, sd in enumerate(seed_dirs):
        cif = sd / "model.cif"
        if not cif.exists():
            continue
        obj_name = f"r{rank:03d}_{short_tpl}_d{int(seq_id):04d}_s{j}"
        cmd.load(str(cif), obj_name)
        cmd.color(color_name, f"{obj_name} and chain A")
        all_design_objs.append(obj_name)
        loaded_by_template[tpl].append(obj_name)
        if (tpl, seq_id) in top10_keys:
            top10_objs.append(obj_name)

    n_loaded = min(len(seed_dirs), N_REPLICATES)
    is_top10 = "TOP10" if (tpl, seq_id) in top10_keys else ""
    print(f"  Rank {rank}: {tpl} design_{seq_id} "
          f"({n_loaded} replicates) {is_top10}")

# ── 3. Create groups ──
if all_design_objs:
    cmd.group("designs", " ".join(all_design_objs))

# Top 10 global group
if top10_objs:
    cmd.group("top10_global", " ".join(top10_objs))

# Per-template design groups
for tpl in templates_seen:
    objs = loaded_by_template.get(tpl, [])
    short = tpl.replace("G4FP_", "")
    if objs:
        cmd.group(f"tpl_{short}", " ".join(objs))

# ── 4. Align everything to first template ──
ref_obj = template_objs[0] if template_objs else (
    all_design_objs[0] if all_design_objs else None)

if ref_obj:
    print("\nAligning structures...")
    all_objs = cmd.get_names("objects")
    for obj in all_objs:
        if obj != ref_obj:
            try:
                cmd.align(obj, ref_obj)
            except Exception:
                pass

# ── 5. Style settings ──
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops", 1)
cmd.set("cartoon_oval_length", 1.2)
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.bg_color("white")

# Base representation: cartoon for protein + DNA
cmd.hide("everything")
cmd.show("cartoon")

# ── 6. DNA (G-quadruplex) as cartoon ──
# Chain B = DNA in AF3 output
cmd.color("gray50", "chain B")

# Create a group selection for all G4 DNA
cmd.select("G4_DNA", "chain B")
cmd.deselect()
print("  G4 DNA (chain B): cartoon")

# ── 7. K+ ions as nb_spheres, colored by atom type ──
# Chain C = K+ ions in AF3 output
cmd.show("nb_spheres", "chain C")
cmd.util.cbag("chain C")  # color by atom type (green for K)
cmd.set("sphere_scale", 0.5, "chain C")

cmd.select("K_ions", "chain C")
cmd.deselect()
print("  K+ ions (chain C): nb_spheres, colored by atom")

# ── 8. Chromophore TYG selection group (all structures) ──
cmd.select("chromophore_TYG", f"({CHROMOPHORE_RESI}) and chain A")
cmd.show("sticks", "chromophore_TYG")
cmd.color("yellow", "chromophore_TYG")
cmd.deselect()
print(f"  Chromophore TYG ({CHROMOPHORE_RESI}): yellow sticks")

# ── 9. Toggle visibility ──
# Templates semi-transparent in background, disabled by default
if all_tpl_objs:
    cmd.disable("templates")

# Zoom to fit
cmd.zoom()
cmd.orient()

# ── 10. Save session ──
output_path = str(base_dir / output_name)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
cmd.save(output_path)

n_total = len(all_design_objs)
print(f"\nPyMOL session saved: {output_path}")
print(f"  {n_total} design objects ({len(selected)} designs x {N_REPLICATES} replicates)")
print(f"  {len(all_tpl_objs)} template objects "
      f"({len(template_objs)} input PDBs + {len(tpl_replicate_objs)} AF3 replicates)")
print(f"\nGroups:")
print(f"  templates      - input PDBs + AF3 replicates (disabled)")
print(f"  designs        - all selected design structures")
print(f"  top10_global   - top 10 overall designs")
for tpl in templates_seen:
    short = tpl.replace("G4FP_", "")
    n = len(loaded_by_template.get(tpl, []))
    print(f"  tpl_{short:20s} - {n} objects")
print(f"\nSelections:")
print(f"  chromophore_TYG - all TYG residues (197-199), yellow sticks")
print(f"  G4_DNA          - G-quadruplex DNA (chain B)")
print(f"  K_ions          - potassium ions (chain C)")
print(f"\nOpen with: pymol {output_name}")

#!/usr/bin/env python3
"""
Build a PyMOL session (.pse) with the selected G4FP design candidates.

Loads all designs from 07_final_candidates.csv, colored by template
using the same tab10 palette as the analysis plots. Chromophore (TYG)
residues are shown as sticks, auto-detected per template.

Usage:
    source /programs/sbgrid.shrc
    pymol -cq build_pymol_session.py
    # or: pymol -cq build_pymol_session.py -- --top 50 --output custom.pse
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict

# Parse args after "--"
top_n = None  # None = load all candidates
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

try:
    _script_path = Path(__file__).resolve()
    base_dir = _script_path.parent
    if not (base_dir / "07_results").exists():
        base_dir = Path.cwd()
except NameError:
    base_dir = Path.cwd()

# ── tab10 colors matching matplotlib plots ──
# Sorted template names → tab10 index, same order as 07_aggregate_and_rank.py
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

# ── Auto-detect chromophore TYG position per template ──


def detect_chromophore_residues(output_dir):
    """Find TYG tripeptide position in the protein sequence.

    Returns 1-indexed (start, end) or None.
    """
    # Strategy 1: LigandMPNN reference FASTA
    fa_dir = output_dir / "01_ligandmpnn" / "seqs"
    if fa_dir.exists():
        for fa in sorted(fa_dir.glob("*.fa")):
            with open(fa) as f:
                lines = f.readlines()
            if len(lines) >= 2:
                seq = lines[1].strip()
                idx = seq.find("TYG")
                if idx >= 0:
                    return (idx + 1, idx + 3)  # 1-indexed
            break

    # Strategy 2: AF3 input JSON
    import json
    for sub in ["02_alphafold3_inputs_bound", "02_alphafold3_inputs_apo"]:
        af3_dir = output_dir / sub
        if not af3_dir.exists():
            continue
        for jf in sorted(af3_dir.glob("*.json")):
            if jf.name.endswith('.bak'):
                continue
            try:
                data = json.load(open(jf))
                seq = data['sequences'][0]['protein']['sequence']
                idx = seq.find("TYG")
                if idx >= 0:
                    return (idx + 1, idx + 3)
            except Exception:
                pass
            break
    return None


# ── Load candidates ──
candidates_csv = base_dir / "07_results" / "07_final_candidates.csv"
if not candidates_csv.exists():
    print(f"ERROR: {candidates_csv} not found. Run 07_aggregate_and_rank.py first.")
    sys.exit(1)

with open(candidates_csv) as f:
    reader = csv.DictReader(f)
    candidates = list(reader)

if top_n is not None:
    candidates = candidates[:top_n]

print(f"\nBuilding PyMOL session with {len(candidates)} designs...")

# ── Discover design structure files ──


def find_best_model(predictions_dir, seq_id, state):
    """Find the best (top-level) model CIF for a design."""
    if not predictions_dir.exists():
        return None

    # Try padded dir format: design_NNNN_state/design_N_state/
    padded = predictions_dir / f"design_{int(seq_id):04d}_{state}"
    if padded.is_dir():
        for sub in padded.iterdir():
            if sub.is_dir():
                cifs = list(sub.glob("*_model.cif"))
                if cifs:
                    return cifs[0]
        cifs = list(padded.glob("*_model.cif"))
        if cifs:
            return cifs[0]

    # Try non-padded: design_N_state/
    non_padded = predictions_dir / f"design_{seq_id}_{state}"
    if non_padded.is_dir():
        cifs = list(non_padded.glob("*_model.cif"))
        if cifs:
            return cifs[0]

    # Try with timestamps
    pattern = re.compile(
        rf'^design_{seq_id}_{re.escape(state)}(?:_\d{{8}}_\d{{6}})?$'
    )
    for d in predictions_dir.iterdir():
        if d.is_dir() and pattern.match(d.name):
            cifs = list(d.glob("*_model.cif"))
            if cifs:
                return cifs[0]
            for sub in d.iterdir():
                if sub.is_dir():
                    cifs = list(sub.glob("*_model.cif"))
                    if cifs:
                        return cifs[0]
    return None


# ── Reinitialize PyMOL ──
cmd.reinitialize()
cmd.set("cif_use_auth", 0)

# ── Detect chromophore per template ──
templates_seen = sorted(set(row["template"] for row in candidates))
template_chromophore = {}  # template_name → "resi X-Z"
for tpl in templates_seen:
    output_dir = base_dir / f"output_{tpl}"
    cr = detect_chromophore_residues(output_dir)
    if cr is not None:
        template_chromophore[tpl] = f"resi {cr[0]}-{cr[1]}"
        print(f"  {tpl}: chromophore TYG at residues {cr[0]}-{cr[1]}")
    else:
        template_chromophore[tpl] = "resi 197-199"
        print(f"  {tpl}: chromophore TYG not found, using 197-199")

# Register template colors in PyMOL
for tpl, rgb in TEMPLATE_COLORS.items():
    cmd.set_color(f"tpl_{tpl}", rgb)

# ── 1. Load template input PDBs ──
template_objs = []
for tpl in templates_seen:
    template_pdb = base_dir / "inputs" / f"{tpl}.pdb"
    if template_pdb.exists():
        obj_name = f"tpl_{tpl}"
        cmd.load(str(template_pdb), obj_name)
        rgb = TEMPLATE_COLORS.get(tpl, [0.7, 0.7, 0.7])
        cmd.set_color(f"c_{obj_name}", rgb)
        cmd.color(f"c_{obj_name}", obj_name)
        cmd.show("cartoon", obj_name)
        cmd.set("cartoon_transparency", 0.6, obj_name)
        template_objs.append(obj_name)
        print(f"  Loaded template: {tpl}")

if template_objs:
    cmd.group("templates", " ".join(template_objs))

# ── 2. Load all candidate designs (holo + apo) ──
loaded_by_template = defaultdict(list)  # template → [obj_names]
loaded_holo = []
loaded_apo = []

for row in candidates:
    seq_id = row["seq_id"]
    tpl = row["template"]
    rank = int(row["final_rank"])
    fitness = float(row["fitness_score"])
    output_dir = base_dir / f"output_{tpl}"
    rgb = TEMPLATE_COLORS.get(tpl, [0.7, 0.7, 0.7])

    for display_state, disk_state in [("holo", "bound"), ("apo", "apo")]:
        pred_dir = output_dir / f"03_alphafold3_predictions_{disk_state}"
        model_cif = find_best_model(pred_dir, seq_id, disk_state)
        if model_cif is None:
            continue

        obj_name = f"r{rank:03d}_{tpl}_d{int(seq_id):04d}_{display_state}"
        cmd.load(str(model_cif), obj_name)

        # Color by template (matching plots)
        cmd.set_color(f"c_{obj_name}", rgb)
        cmd.color(f"c_{obj_name}", obj_name)
        cmd.show("cartoon", obj_name)

        if display_state == "holo":
            loaded_holo.append(obj_name)
        else:
            loaded_apo.append(obj_name)
        loaded_by_template[tpl].append(obj_name)

    print(f"  Rank {rank}: {tpl} design_{seq_id} (fitness={fitness:.3f})")

# ── 3. Create groups ──
if loaded_holo:
    cmd.group("designs_holo", " ".join(loaded_holo))
if loaded_apo:
    cmd.group("designs_apo", " ".join(loaded_apo))

# Per-template groups
for tpl in templates_seen:
    objs = loaded_by_template.get(tpl, [])
    if objs:
        cmd.group(f"tpl_grp_{tpl}", " ".join(objs))

# ── 4. Align everything to the first template ──
ref_obj = template_objs[0] if template_objs else (
    loaded_holo[0] if loaded_holo else None)

if ref_obj:
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

# Show chain A cartoon only
cmd.hide("everything")
cmd.show("cartoon", "chain A")

# ── 6. Show chromophore TYG as sticks per template ──
for tpl in templates_seen:
    resi_sel = template_chromophore[tpl]
    objs = loaded_by_template.get(tpl, [])
    # Include the template input PDB object
    tpl_obj = f"tpl_{tpl}"
    if tpl_obj in cmd.get_names("objects"):
        objs = [tpl_obj] + objs

    for obj in objs:
        sel_name = f"chrom_{obj}"
        cmd.select(sel_name, f"{obj} and {resi_sel} and chain A")
        cmd.show("sticks", sel_name)
        cmd.color("yellow", sel_name)
        cmd.deselect()

# ── 7. Toggle visibility ──
# Start with apo hidden
if loaded_apo:
    cmd.disable("designs_apo")

# Templates semi-transparent in background
if template_objs:
    cmd.disable("templates")

# Zoom to fit
cmd.zoom("chain A")
cmd.orient()

# ── 8. Save session ──
output_path = str(base_dir / output_name)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
cmd.save(output_path)

n_total = len(loaded_holo) + len(loaded_apo)
print(f"\nPyMOL session saved: {output_path}")
print(f"  {len(loaded_holo)} holo + {len(loaded_apo)} apo structures")
print(f"  {len(templates_seen)} templates (colored by tab10 palette)")
print(f"  Groups: templates, designs_holo, designs_apo, tpl_grp_<name>")
print(f"  Chromophore TYG shown as yellow sticks (per-template positions)")
print(f"  Apo + templates disabled by default (toggle in PyMOL)")
print(f"\nOpen with: pymol {output_name}")

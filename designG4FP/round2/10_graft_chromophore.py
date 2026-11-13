#!/usr/bin/env python3
"""
10_graft_chromophore.py -- rebuild the mature chromophore into AF3 models and
convert them to PDB for LigandMPNN.

Why
---
AF3 does not model the mature GFP chromophore.  It predicts the Thr-Tyr-Gly
triad as three ordinary amino acids in an extended conformation: across the 14
round-2 input structures the Thr-Ca to Gly-Ca distance is 4.66-5.35 A, versus
3.10 A in the cyclized CRO of the design templates.  Any metric computed on
that triad therefore describes an unmatured loop, not a chromophore.

This is the same failure mode as the Zn(II)-only models in Chen, Wu et al.
(Nature 2026, doi 10.1038/s41586-026-11167-2): with the mechanistically
essential species absent, active and inactive designs were indistinguishable
(3/89 hit rate); once the catalytic water was modelled explicitly the hit rate
went to 31/38.  For G4FP the mechanistically essential species is the mature
chromophore -- fluorescence IS the chromophore.

What it does
------------
For each design model, superpose its parent template (which carries CRO as a
HETATM group) onto the model using Ca atoms of shared residues within a radius
of the chromophore, then transplant the CRO group in place of the model's TYG
triad.  Local superposition rather than whole-chain, because the designed
termini differ by ~11 A RMSD while the barrel around the chromophore agrees to
~1.3 A.

Output is PDB, not mmCIF: LigandMPNN's run.py reads PDB only (ProDy
parsePDB), so this doubles as the format conversion.

Usage
-----
    python 10_graft_chromophore.py --models  structure_files/designs \\
                                   --templates ../inputs \\
                                   --out inputs_round2 --radius 15

Diagnostics are written to <out>/10_graft_report.csv.  Inspect
`graft_fit_rmsd` and `new_clashes` before using any output: a design whose
barrel cannot accommodate the cyclized chromophore is itself informative and
should not be silently designed on.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}
CHROMOPHORE_MOTIF = re.compile(r"VTTL(TYG)VQ")
CHROMOPHORE_MOTIF_LOOSE = re.compile(r"TTL(TYG)V")
DNA_COMPS = {"DA", "DG", "DC", "DT"}
# CRO atom names whose positions correspond to the Ca of Thr, Tyr, Gly
CRO_ANCHORS = ("CA1", "CA2", "CA3")


# ---------------------------------------------------------------- parsing

def parse_cif(path: Path) -> List[dict]:
    cols: List[str] = []
    out: List[dict] = []
    for line in path.read_text(errors="replace").splitlines():
        s = line.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1])
            continue
        if not cols:
            continue
        if s.startswith(("ATOM", "HETATM")):
            f = s.split()
            if len(f) < len(cols):
                continue
            d = dict(zip(cols, f))
            sid = d.get("label_seq_id", ".")
            if sid in (".", "?"):
                sid = d.get("auth_seq_id", ".")
            if sid in (".", "?"):
                continue
            try:
                out.append({"record": s.split()[0], "chain": d["label_asym_id"],
                            "resi": int(sid), "resn": d["label_comp_id"],
                            "name": d["label_atom_id"].strip('"'),
                            "element": d.get("type_symbol", ""),
                            "xyz": np.array([float(d["Cartn_x"]), float(d["Cartn_y"]),
                                             float(d["Cartn_z"])]),
                            "b": float(d.get("B_iso_or_equiv", 0.0))})
            except (KeyError, ValueError):
                continue
        elif out and s.startswith("#"):
            break
    return out


def parse_pdb(path: Path) -> List[dict]:
    out = []
    for l in path.read_text(errors="replace").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            try:
                out.append({"record": l[:6].strip(), "chain": l[21],
                            "resi": int(re.sub(r"\D", "", l[22:27]) or 0),
                            "resn": l[17:20].strip(), "name": l[12:16].strip(),
                            "element": l[76:78].strip(),
                            "xyz": np.array([float(l[30:38]), float(l[38:46]),
                                             float(l[46:54])]),
                            "b": float(l[60:66] or 0)})
            except ValueError:
                continue
    return out


def load_any(path: Path) -> List[dict]:
    return parse_cif(path) if path.suffix.lower() == ".cif" else parse_pdb(path)


# ---------------------------------------------------------------- geometry

def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t mapping P onto Q (proper rotation only)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Qc - R @ Pc


def protein_chain(atoms: List[dict]) -> str:
    counts: Dict[str, set] = {}
    for a in atoms:
        if a["resn"] in THREE2ONE:
            counts.setdefault(a["chain"], set()).add(a["resi"])
    return max(counts, key=lambda c: len(counts[c]))


def find_triad(atoms: List[dict], chain: str) -> Optional[List[int]]:
    comp = {a["resi"]: a["resn"] for a in atoms if a["chain"] == chain}
    order = sorted(comp)
    seq = "".join(THREE2ONE.get(comp[s], "X") for s in order)
    for pat in (CHROMOPHORE_MOTIF, CHROMOPHORE_MOTIF_LOOSE):
        m = pat.search(seq)
        if m:
            i = m.start(1)
            return [order[i], order[i + 1], order[i + 2]]
    return None


# ---------------------------------------------------------------- grafting

def graft(model_atoms: List[dict], template_atoms: List[dict],
          radius: float = 15.0) -> Tuple[List[dict], dict]:
    mchain = protein_chain(model_atoms)
    tchain = protein_chain(template_atoms)
    triad = find_triad(model_atoms, mchain)
    if triad is None:
        raise ValueError("chromophore triad not found in model")

    cro = {a["name"]: a for a in template_atoms
           if a["chain"] == tchain and a["resn"] == "CRO"}
    if not all(k in cro for k in CRO_ANCHORS):
        raise ValueError("template has no CRO group with CA1/CA2/CA3")
    cro_centre = np.mean([cro[k]["xyz"] for k in CRO_ANCHORS], axis=0)

    tca = {a["resi"]: a["xyz"] for a in template_atoms
           if a["chain"] == tchain and a["name"] == "CA"}
    mca = {a["resi"]: a["xyz"] for a in model_atoms
           if a["chain"] == mchain and a["name"] == "CA"}
    # The triad and the residue immediately after it have no 1:1 counterpart
    # (the template collapses TYG into a single CRO residue).
    excl = set(triad) | {triad[-1] + 1}
    shared = [n for n in sorted(set(tca) & set(mca) - excl)
              if np.linalg.norm(tca[n] - cro_centre) < radius]
    if len(shared) < 20:
        raise ValueError(f"only {len(shared)} shared residues within {radius} A")

    P = np.array([tca[n] for n in shared])
    Q = np.array([mca[n] for n in shared])
    R, t = kabsch(P, Q)
    fit_rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean()))

    # keep everything except the model's triad residues
    kept = [a for a in model_atoms
            if not (a["chain"] == mchain and a["resi"] in triad)]
    cro_atoms = []
    for name, a in cro.items():
        cro_atoms.append({"record": "HETATM", "chain": mchain, "resi": triad[1],
                          "resn": "CRO", "name": name,
                          "element": a["element"] or name[0],
                          "xyz": (R @ a["xyz"]) + t, "b": 1.00})

    # clash check: grafted CRO against everything kept, excluding its own
    # covalent neighbours (residues flanking the triad)
    flank = {triad[0] - 1, triad[-1] + 1}
    others = np.array([a["xyz"] for a in kept
                       if not (a["chain"] == mchain and a["resi"] in flank)])
    C = np.array([a["xyz"] for a in cro_atoms])
    d = np.sqrt(((C[:, None, :] - others[None, :, :]) ** 2).sum(-1))
    report = {"triad": triad, "n_superposed": len(shared),
              "graft_fit_rmsd": round(fit_rmsd, 3),
              "cro_min_dist": round(float(d.min()), 2),
              "new_clashes": int((d < 2.2).sum())}
    return kept + cro_atoms, report


# ---------------------------------------------------------------- output

def write_pdb(atoms: List[dict], path: Path) -> None:
    order = {"ATOM": 0, "HETATM": 1}
    atoms = sorted(atoms, key=lambda a: (a["chain"], order.get(a["record"], 1), a["resi"]))
    with open(path, "w") as fh:
        for i, a in enumerate(atoms, 1):
            name = a["name"]
            nm = f" {name:<3}" if len(name) < 4 else name
            x, y, z = a["xyz"]
            fh.write(f"{a['record']:<6}{i:>5} {nm}{'':1}{a['resn']:>3} {a['chain']}"
                     f"{a['resi']:>4}{'':4}{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                     f"{1.00:>6.2f}{a['b']:>6.2f}{'':10}{a['element']:>2}\n")
        fh.write("END\n")


FLAT_RE = re.compile(r"_d(\d{3,5})_s(\d+)\.(cif|pdb)$")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", type=Path, required=True,
                    help="directory tree of design models (.cif or .pdb)")
    ap.add_argument("--templates", type=Path, required=True,
                    help="directory of parent template PDBs carrying CRO")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radius", type=float, default=15.0,
                    help="superposition radius around the chromophore (A)")
    ap.add_argument("--seed-only", default="s0",
                    help="model suffix to convert, or 'all'")
    args = ap.parse_args()

    templates = {}
    for p in sorted(args.templates.glob("*.pdb")):
        templates[p.stem.replace("G4FP_", "")] = p
    if not templates:
        sys.exit(f"no template PDBs in {args.templates}")
    print(f"{len(templates)} templates loaded")

    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    models = [p for p in args.models.rglob("*") if p.suffix.lower() in (".cif", ".pdb")]
    for m in sorted(models):
        mm = FLAT_RE.search(m.name)
        if not mm:
            continue
        if args.seed_only != "all" and not m.name.endswith(f"_{args.seed_only}{m.suffix}"):
            continue
        tpl_key = m.parent.name
        if tpl_key not in templates:
            print(f"  SKIP {m.name}: no template named {tpl_key}")
            continue
        try:
            atoms, rep = graft(load_any(m), load_any(templates[tpl_key]), args.radius)
        except ValueError as e:
            print(f"  FAIL {m.name}: {e}")
            rows.append({"model": m.name, "template": tpl_key, "error": str(e)})
            continue
        out = args.out / f"{tpl_key}_d{mm.group(1)}_CRO.pdb"
        write_pdb(atoms, out)
        rep.update({"model": m.name, "template": tpl_key, "output": out.name})
        rows.append(rep)
        print(f"  {out.name}: fit {rep['graft_fit_rmsd']} A over {rep['n_superposed']} res, "
              f"CRO min dist {rep['cro_min_dist']} A, clashes {rep['new_clashes']}")

    df = pd.DataFrame(rows)
    df.to_csv(args.out / "10_graft_report.csv", index=False)
    ok = df[df.get("new_clashes", pd.Series(dtype=float)).notna()] if len(df) else df
    print(f"\nWrote {len(ok)} grafted structures to {args.out}/")
    if len(ok):
        print(f"graft_fit_rmsd: {ok.graft_fit_rmsd.min():.2f}-{ok.graft_fit_rmsd.max():.2f} A")
        print(f"new clashes (<2.2 A): {ok.new_clashes.min()}-{ok.new_clashes.max()}")


if __name__ == "__main__":
    main()

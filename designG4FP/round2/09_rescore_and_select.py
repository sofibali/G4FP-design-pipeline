#!/usr/bin/env python3
"""
09_rescore_and_select.py -- G4FP round-2 rescoring and input selection.

Why this script exists
----------------------
Three defects were found in the round-1 scoring (see round2/README.md):

  1. `compute_per_residue_plddt()` in 05_analyze_af3_structures.py maps
     atom-level pLDDT onto residues by assuming a UNIFORM atoms-per-residue
     (len(chain_plddts) / n_residues).  Real residues have 4-14 atoms, so the
     mapping drifts cumulatively along the chain.  Combined with the hardcoded
     chromophore range 197-199, the reported chromophore pLDDT for the
     r28_l2_cro_mod2 lineage actually averaged residues ~201-204
     (Cys-Phe-Ser-Arg).  Measured error on 14 structures: +12.8 pLDDT points
     on average, range 6.6-19.8.

  2. The chromophore triad is not at 197-199 in most constructs.  Across the
     14 tested-variant structures it starts at residue 187, 195, 196 or 197.

  3. `apo_plddt_vs_template` was NaN for every design in the committed
     results, so score_component_rel_apo was a constant 0.5 and the apo term
     did no ranking work at all.

This script recomputes the chromophore metrics correctly, adds a chromophore
geometry gate, recomputes the fitness score, takes the Pareto frontier and
selects the top N designs per starting template as LigandMPNN round-2 inputs.

Usage
-----
    # Sofia's AF3 output tree (seed-*_sample-* subdirectories)
    python 09_rescore_and_select.py --af3-root /path/to/predictions_bound \\
        --round1-table ../07_results/plot_data_csvs/07_all_designs_scored.csv \\
        --out 09_results --top-n 5

    # Flat CIF export (one file per seed, *_d<NNNN>_s<k>.cif)
    python 09_rescore_and_select.py --cif-root structure_files/designs \\
        --layout flat --out 09_results --top-n 5

    # Check how the experimentally tested variants rank under the weights
    python 09_rescore_and_select.py ... --validate tested_variants.csv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}

# The chromophore-forming triad sits inside this invariant cpEGFP stretch in
# every G4FP construct.  Matching on sequence rather than residue number makes
# the script independent of construct length and linker choice.
CHROMOPHORE_MOTIF = re.compile(r"VTTL(TYG)VQ")
CHROMOPHORE_MOTIF_LOOSE = re.compile(r"TTL(TYG)V")

DNA_COMPS = {"DA", "DG", "DC", "DT"}


# --------------------------------------------------------------------------
# mmCIF parsing
# --------------------------------------------------------------------------

def parse_cif_atoms(path: Path) -> List[tuple]:
    """Minimal mmCIF atom_site reader.

    Returns (chain, seq_id, comp_id, atom_id, x, y, z, bfactor) tuples.
    AlphaFold3 writes per-atom pLDDT into B_iso_or_equiv.
    """
    cols: List[str] = []
    atoms: List[tuple] = []
    with open(path) as fh:
        for line in fh:
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
                    atoms.append((
                        d["label_asym_id"], int(sid), d["label_comp_id"],
                        d["label_atom_id"],
                        float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"]),
                        float(d.get("B_iso_or_equiv", "nan")),
                    ))
                except (KeyError, ValueError):
                    continue
            elif atoms and s.startswith("#"):
                break
    return atoms


def protein_chain_id(atoms: List[tuple]) -> Optional[str]:
    """The protein chain is the one with the most amino-acid residues."""
    counts: Dict[str, set] = defaultdict(set)
    for ch, sid, comp, *_ in atoms:
        if comp in THREE2ONE:
            counts[ch].add(sid)
    if not counts:
        return None
    return max(counts, key=lambda c: len(counts[c]))


# --------------------------------------------------------------------------
# Per-residue pLDDT -- the corrected version
# --------------------------------------------------------------------------

def per_residue_plddt(atoms: List[tuple], chain: str) -> Dict[int, float]:
    """Mean pLDDT per residue, grouped by ACTUAL residue id.

    This is the fix for defect (1).  The original implementation binned the
    atom array into equal-sized slices, which silently shifts the mapping by
    several residues once the chain is a few hundred residues long.
    """
    acc: Dict[int, List[float]] = defaultdict(list)
    for ch, sid, comp, aid, x, y, z, b in atoms:
        if ch == chain and not np.isnan(b):
            acc[sid].append(b)
    return {sid: float(np.mean(v)) for sid, v in acc.items()}


def chain_sequence(atoms: List[tuple], chain: str) -> Tuple[str, List[int]]:
    comp: Dict[int, str] = {}
    for ch, sid, cp, *_ in atoms:
        if ch == chain:
            comp[sid] = cp
    order = sorted(comp)
    seq = "".join(THREE2ONE.get(comp[s], "X") for s in order)
    return seq, order


def find_chromophore_triad(atoms: List[tuple], chain: str,
                           fallback_range: Optional[Tuple[int, int]] = None
                           ) -> Optional[List[int]]:
    """Locate the Thr-Tyr-Gly chromophore triad; returns three residue ids.

    Sequence motif first (defect 2), then an explicit CRO heteroresidue if the
    structure carries a mature chromophore, then an optional hardcoded range.
    """
    seq, order = chain_sequence(atoms, chain)
    for pattern in (CHROMOPHORE_MOTIF, CHROMOPHORE_MOTIF_LOOSE):
        m = pattern.search(seq)
        if m:
            i = m.start(1)
            return [order[i], order[i + 1], order[i + 2]]
    cro = sorted({sid for ch, sid, comp, *_ in atoms if ch == chain and comp == "CRO"})
    if cro:
        return [cro[0], cro[0], cro[0]]
    if fallback_range:
        lo, hi = fallback_range
        present = [s for s in order if lo <= s <= hi]
        if len(present) >= 3:
            return present[:3]
    return None


def _atom_xyz(atoms: List[tuple], chain: str, sid: int, name: str) -> Optional[np.ndarray]:
    for ch, s, comp, aid, x, y, z, b in atoms:
        if ch == chain and s == sid and aid == name:
            return np.array([x, y, z])
    return None


def cyclization_distance(atoms: List[tuple], chain: str, triad: List[int]) -> float:
    """Thr carbonyl C to Gly amide N -- the bond formed on chromophore maturation.

    In a mature CRO this is a covalent bond (~1.3 A).  In AF3 models of the
    unmatured TYG triad it sits near 3 A; a grossly distorted triad reads much
    longer.  Used as a geometry GATE, not as a ranking term: across the 14
    round-2 input structures it spans only 2.90-3.20 A, far too narrow to rank
    on.
    """
    if len(set(triad)) == 1:          # mature CRO modelled as a single residue
        return float("nan")
    c = _atom_xyz(atoms, chain, triad[0], "C")
    n = _atom_xyz(atoms, chain, triad[2], "N")
    if c is None or n is None:
        return float("nan")
    return float(np.linalg.norm(c - n))


# --------------------------------------------------------------------------
# Per-model metric extraction
# --------------------------------------------------------------------------

def metrics_for_model(cif: Path, fallback_range: Optional[Tuple[int, int]]) -> Optional[dict]:
    atoms = parse_cif_atoms(cif)
    if not atoms:
        return None
    chain = protein_chain_id(atoms)
    if chain is None:
        return None
    plddt = per_residue_plddt(atoms, chain)
    if not plddt:
        return None
    triad = find_chromophore_triad(atoms, chain, fallback_range)
    rec = {
        "mean_plddt": float(np.mean(list(plddt.values()))),
        "min_plddt": float(np.min(list(plddt.values()))),
        "n_residues": len(plddt),
        "chromophore_plddt": float("nan"),
        "chromophore_start": np.nan,
        "d_cyclization": float("nan"),
    }
    if triad:
        uniq = list(dict.fromkeys(triad))
        rec["chromophore_plddt"] = float(np.mean([plddt[t] for t in uniq]))
        rec["chromophore_start"] = int(triad[0])
        rec["d_cyclization"] = cyclization_distance(atoms, chain, triad)
    # protein-DNA closest approach: catches physically invalid interfaces
    dna = np.array([a[4:7] for a in atoms if a[2] in DNA_COMPS])
    if len(dna):
        prot = np.array([a[4:7] for a in atoms if a[0] == chain])
        dmat = np.sqrt(((prot[:, None, :] - dna[None, :, :]) ** 2).sum(-1))
        rec["min_protein_dna_dist"] = float(dmat.min())
        rec["n_clash_pairs"] = int((dmat < 2.5).sum())
    return rec


def read_summary_confidences(d: Path) -> dict:
    """Pick up pTM / ipTM / ranking_score from whatever AF3 wrote next to the model."""
    out: dict = {}
    for f in sorted(d.glob("*summary_confidences*.json")):
        try:
            c = json.load(open(f))
        except Exception:
            continue
        for k in ("ptm", "iptm", "ranking_score", "fraction_disordered"):
            if k in c and k not in out:
                out[k] = c[k]
        break
    return out


# --------------------------------------------------------------------------
# Input discovery
# --------------------------------------------------------------------------

FLAT_RE = re.compile(r"_d(\d{3,5})_s(\d+)\.cif$")


def norm_id(gid: str) -> str:
    """Normalise a design identifier for joining across sources.

    Flat CIF exports name the directory `r28_l2_cro_mod2` while the round-1
    tables use `G4FP_r28_l2_cro_mod2_design_0072`, and design indices appear
    both zero-padded and bare.  Strip the prefix and pad the index.
    """
    g = str(gid)
    if g.startswith("G4FP_"):
        g = g[len("G4FP_"):]
    m = re.search(r"^(.*)_design_(\d+)$", g)
    if m:
        g = f"{m.group(1)}_design_{int(m.group(2)):04d}"
    return g


def discover(root: Path, layout: str) -> Dict[Tuple[str, str], List[Path]]:
    """Return {(template, design_id): [model cif, ...]}."""
    found: Dict[Tuple[str, str], List[Path]] = defaultdict(list)
    if layout in ("flat", "auto"):
        for cif in root.rglob("*.cif"):
            m = FLAT_RE.search(cif.name)
            if m:
                found[(cif.parent.name, m.group(1))].append(cif)
        if found and layout == "auto":
            return found
    if layout in ("af3", "auto"):
        for d in sorted(p for p in root.rglob("*") if p.is_dir()):
            if not re.match(r"^seed-\d+_sample-\d+$", d.name):
                continue
            models = sorted(d.glob("*model*.cif"))
            if not models:
                continue
            design_dir = d.parent
            m = re.search(r"design_(\d+)", design_dir.name)
            design_id = m.group(1) if m else design_dir.name
            template = design_dir.parent.name.replace("output_", "")
            found[(template, design_id)].extend(models)
    return found


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

DEFAULT_WEIGHTS = {
    # Interface quality: the sensor must actually engage the G4.
    "iptm": 0.30,
    # Chromophore order in the holo state -- now measured on the real triad.
    "chrom_holo": 0.25,
    # Global fold of the bound state.
    "holo_ptm": 0.15,
    # Apo order relative to the parent template.  Round 1 weighted this so that
    # LOWER was better (the switch hypothesis).  The only two experimentally
    # fluorescent variants sit ABOVE their template (+0.65, +0.69) while all 19
    # dark ones sit below (-2.66 to -0.005), so the sign is inverted here.
    # n = 2, and the metric is partly a property of the template, so this is a
    # hypothesis encoded for testing, not an established effect.  Use
    # --rel-apo-sign -1 to restore the round-1 direction, or --w-rel-apo 0.
    "rel_apo": 0.15,
    # Holo/apo conformational change.
    "rmsd": 0.10,
    # Agreement across seeds.
    "confidence": 0.05,
}


def minmax(s: pd.Series) -> pd.Series:
    lo, hi = s.min(), s.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        return pd.Series(0.5, index=s.index)
    return (s - lo) / (hi - lo)


def compute_fitness(df: pd.DataFrame, weights: dict, rel_apo_sign: float = 1.0) -> pd.DataFrame:
    df = df.copy()
    nan = pd.Series(np.nan, index=df.index)
    comp = {
        "iptm": minmax(df.get("holo_iptm", nan)),
        "chrom_holo": minmax(df.get("holo_chromophore_plddt_fixed", nan)),
        "holo_ptm": minmax(df.get("holo_ptm", nan)),
        "rmsd": minmax(df.get("global_rmsd", nan)),
        "confidence": 1.0 - minmax(df.get("holo_chromophore_plddt_fixed_sd", nan)),
    }
    if "apo_plddt_vs_template" in df and df["apo_plddt_vs_template"].notna().any():
        comp["rel_apo"] = minmax(rel_apo_sign * df["apo_plddt_vs_template"])
    else:
        comp["rel_apo"] = pd.Series(0.5, index=df.index)
        print("  NOTE: apo_plddt_vs_template unavailable -- the rel_apo term is "
              "inert. This was silently true throughout round 1.")

    total = 0.0
    for k, w in weights.items():
        col = f"score_component_{k}"
        df[col] = comp[k].fillna(0.5)
        total = total + w * df[col]
    df["fitness_score"] = total
    return df


def apply_gates(df: pd.DataFrame, min_holo_plddt: float, min_holo_ptm: float,
                min_holo_iptm: float, cyc_min: float, cyc_max: float,
                max_clash_pairs: int) -> pd.DataFrame:
    n0 = len(df)
    mask = pd.Series(True, index=df.index)
    reasons: Dict[str, int] = {}

    def gate(name, m):
        nonlocal mask
        reasons[name] = int((~m & mask).sum())
        mask = mask & m

    if "holo_mean_plddt" in df:
        gate(f"holo_mean_plddt >= {min_holo_plddt}", df["holo_mean_plddt"] >= min_holo_plddt)
    if "holo_ptm" in df:
        gate(f"holo_ptm >= {min_holo_ptm}", df["holo_ptm"] >= min_holo_ptm)
    if "holo_iptm" in df:
        gate(f"holo_iptm >= {min_holo_iptm}", df["holo_iptm"] >= min_holo_iptm)
    if "d_cyclization" in df:
        d = df["d_cyclization"]
        gate(f"chromophore geometry {cyc_min}-{cyc_max} A",
             ~((d < cyc_min) | (d > cyc_max)) | d.isna())
    # NOTE: every round-1 AF3 holo model carries protein-DNA steric clashes
    # (10-37 atom pairs below 2.5 A per model; templates included, minimum
    # separations 0.7-1.8 A).  A zero-clash gate therefore removes the entire
    # dataset, so this gate is OFF unless a threshold is passed explicitly.
    if "n_clash_pairs" in df and max_clash_pairs is not None:
        gate(f"protein-DNA clash pairs <= {max_clash_pairs}",
             df["n_clash_pairs"].fillna(0) <= max_clash_pairs)

    print(f"\nHard gates: {n0} -> {int(mask.sum())} designs")
    for k, v in reasons.items():
        print(f"    {k:44} removed {v}")
    return df[mask].copy()


def pareto_frontier(df: pd.DataFrame, objectives: List[str]) -> pd.Series:
    """Boolean mask of non-dominated rows (all objectives maximised)."""
    cols = [c for c in objectives if c in df.columns]
    if not cols:
        return pd.Series(True, index=df.index)
    V = df[cols].astype(float)
    V = V.fillna(V.min()).to_numpy()
    n = len(V)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        dominated = np.all(V >= V[i], axis=1) & np.any(V > V[i], axis=1)
        if dominated.any():
            keep[i] = False
    return pd.Series(keep, index=df.index)


def select_per_template(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Top n designs per starting template; Pareto members first, then by fitness."""
    picks = []
    for tpl, g in df.groupby("template"):
        g = g.copy()
        if "on_pareto" not in g:
            g["on_pareto"] = False
        g = g.sort_values(["on_pareto", "fitness_score"], ascending=[False, False])
        picks.append(g.head(n))
    return pd.concat(picks, ignore_index=True) if picks else df.head(0)


# --------------------------------------------------------------------------
# Collection
# --------------------------------------------------------------------------

def collect(root: Path, layout: str, state: str,
            fallback_range: Optional[Tuple[int, int]]) -> pd.DataFrame:
    groups = discover(root, layout)
    if not groups:
        sys.exit(f"No structures found under {root} (layout={layout}).")
    print(f"Found {len(groups)} designs with structures under {root}")
    rows = []
    for (template, design_id), cifs in sorted(groups.items()):
        per_model = []
        for cif in cifs:
            m = metrics_for_model(cif, fallback_range)
            if m:
                m.update(read_summary_confidences(cif.parent))
                per_model.append(m)
        if not per_model:
            continue
        M = pd.DataFrame(per_model)
        rec = {"template": template, "design_id": design_id,
               "global_id": f"{template}_design_{design_id}",
               "n_models": len(M)}
        mapping = [("mean_plddt", f"{state}_mean_plddt"),
                   ("chromophore_plddt", f"{state}_chromophore_plddt_fixed"),
                   ("d_cyclization", "d_cyclization"),
                   ("ptm", f"{state}_ptm"), ("iptm", f"{state}_iptm"),
                   ("min_protein_dna_dist", "min_protein_dna_dist"),
                   ("n_clash_pairs", "n_clash_pairs")]
        for col, out in mapping:
            if col in M:
                rec[out] = float(M[col].mean())
                rec[f"{out}_sd"] = float(M[col].std(ddof=0))
        if "chromophore_start" in M and M["chromophore_start"].notna().any():
            rec["chromophore_start"] = int(M["chromophore_start"].dropna().mode().iloc[0])
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--af3-root", type=Path, help="AF3 prediction tree (seed-*_sample-*)")
    src.add_argument("--cif-root", type=Path, help="Directory tree of exported model CIFs")
    ap.add_argument("--layout", choices=["auto", "af3", "flat"], default="auto")
    ap.add_argument("--apo-root", type=Path, default=None,
                    help="Optional matching apo tree; enables apo and delta metrics")
    ap.add_argument("--round1-table", type=Path, default=None,
                    help="07_all_designs_scored.csv -- merged for RMSD / apo_plddt_vs_template")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=5, help="designs selected per template")
    ap.add_argument("--chromophore-fallback", default=None,
                    help="e.g. 197,199 -- used only if the sequence motif fails")
    ap.add_argument("--min-holo-plddt", type=float, default=70.0)
    ap.add_argument("--min-holo-ptm", type=float, default=0.5)
    ap.add_argument("--min-holo-iptm", type=float, default=0.3)
    ap.add_argument("--cyc-min", type=float, default=2.5,
                    help="minimum Thr-C to Gly-N distance (A); below this the triad "
                         "is sterically impossible")
    ap.add_argument("--cyc-max", type=float, default=3.5,
                    help="maximum Thr-C to Gly-N distance (A); above this the triad "
                         "is too distorted to cyclise")
    ap.add_argument("--max-clash-pairs", type=int, default=None,
                    help="max protein-DNA atom pairs below 2.5 A (default: gate off; "
                         "every round-1 AF3 model has 10-37, so 0 removes everything)")
    for k, v in DEFAULT_WEIGHTS.items():
        ap.add_argument(f"--w-{k.replace('_', '-')}", type=float, default=v)
    ap.add_argument("--rel-apo-sign", type=float, default=1.0,
                    help="+1 rewards apo MORE ordered than template (round-1 data); "
                         "-1 restores the round-1 switch assumption")
    ap.add_argument("--validate", type=Path, default=None,
                    help="CSV with global_id,outcome -- reports ranks of tested variants")
    args = ap.parse_args()

    root = args.af3_root or args.cif_root
    fallback = None
    if args.chromophore_fallback:
        lo, hi = map(int, args.chromophore_fallback.split(","))
        fallback = (lo, hi)

    df = collect(root, args.layout, "holo", fallback)
    print(f"Collected metrics for {len(df)} designs across "
          f"{df['template'].nunique()} templates")

    if args.apo_root:
        apo = collect(args.apo_root, args.layout, "apo", fallback)
        apo = apo.drop(columns=[c for c in ("template", "n_models") if c in apo.columns])
        df = df.merge(apo, on=["global_id", "design_id"], how="left", suffixes=("", "_apo"))
        if "apo_chromophore_plddt_fixed" in df:
            df["chromophore_plddt_diff_fixed"] = (df["holo_chromophore_plddt_fixed"]
                                                  - df["apo_chromophore_plddt_fixed"])

    if args.round1_table and args.round1_table.exists():
        r1 = pd.read_csv(args.round1_table)
        bring = [c for c in ("global_id", "global_rmsd", "apo_plddt_vs_template",
                             "apo_mean_plddt", "sequence",
                             # flat CIF exports carry no summary_confidences JSON,
                             # so pTM / ipTM have to come from the round-1 table
                             "holo_ptm", "holo_iptm", "holo_mean_plddt")
                 if c in r1.columns]
        r1 = r1[bring].copy()
        # Directory names drop the "G4FP_" prefix that the round-1 table carries.
        r1["_key"] = r1["global_id"].map(norm_id)
        df["_key"] = df["global_id"].map(norm_id)
        dup = [c for c in r1.columns if c in df.columns and c not in ("global_id", "_key")]
        if dup:
            print(f"  (structure-derived values kept for: {', '.join(dup)})")
            r1 = r1.drop(columns=dup)
        df = df.merge(r1.drop(columns="global_id"), on="_key", how="left").drop(columns="_key")
        n_matched = int(df[[c for c in bring if c != "global_id"][0]].notna().sum()) \
            if len(bring) > 1 else 0
        print(f"Merged {len(bring) - 1} columns from {args.round1_table.name} "
              f"({n_matched}/{len(df)} designs matched)")

    weights = {k: getattr(args, f"w_{k}") for k in DEFAULT_WEIGHTS}
    print("\nWeights:", {k: round(v, 3) for k, v in weights.items()},
          f"| rel_apo_sign={args.rel_apo_sign:+.0f}")
    df = compute_fitness(df, weights, args.rel_apo_sign)

    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "09_all_designs_rescored.csv", index=False)

    gated = apply_gates(df, args.min_holo_plddt, args.min_holo_ptm, args.min_holo_iptm,
                        args.cyc_min, args.cyc_max, args.max_clash_pairs)
    if gated.empty:
        sys.exit("All designs removed by the hard gates -- loosen them and rerun.")

    objectives = [c for c in ("holo_iptm", "holo_chromophore_plddt_fixed",
                              "holo_ptm", "global_rmsd") if c in gated.columns]
    gated["on_pareto"] = pareto_frontier(gated, objectives)
    print(f"Pareto frontier over {objectives}: {int(gated['on_pareto'].sum())} designs")
    gated.to_csv(args.out / "09_gated_with_pareto.csv", index=False)

    sel = select_per_template(gated, args.top_n)
    sel = sel.sort_values(["template", "fitness_score"], ascending=[True, False])
    sel.to_csv(args.out / "09_round2_ligandmpnn_inputs.csv", index=False)
    print(f"\nSelected {len(sel)} designs "
          f"(up to {args.top_n} per template x {sel['template'].nunique()} templates)")
    print(sel.groupby("template").size().to_string())

    if args.validate and args.validate.exists():
        v = pd.read_csv(args.validate)
        v["_key"] = v["global_id"].map(norm_id)
        df["_key"] = df["global_id"].map(norm_id)
        merged = df.merge(v.drop(columns="global_id"), on="_key", how="inner")
        if merged.empty:
            print("\nValidation file shares no global_id with the rescored set.")
        else:
            merged = merged.sort_values("fitness_score", ascending=False).reset_index(drop=True)
            merged["rank"] = merged.index + 1
            print(f"\nValidation on {len(merged)} experimentally tested variants:")
            cols = [c for c in ("rank", "global_id", "outcome", "fitness_score",
                                "holo_chromophore_plddt_fixed", "d_cyclization")
                    if c in merged.columns]
            print(merged[cols].to_string(index=False))
            merged.to_csv(args.out / "09_validation_tested_variants.csv", index=False)

    print(f"\nWrote results to {args.out}/")


if __name__ == "__main__":
    main()

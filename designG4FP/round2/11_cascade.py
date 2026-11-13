#!/usr/bin/env python3
"""
11_cascade.py -- staged gate cascade for G4FP design selection.

Supersedes the weighted composite in 09_rescore_and_select.py.

Rationale
---------
Modelled on the de novo metalloprotease pipeline of Chen, Wu et al.
(Nature 2026, doi 10.1038/s41586-026-11167-2, Supplementary Methods
"Filtering").  Three properties of that pipeline are adopted here:

  1. NON-COMPENSATORY GATES, NOT A WEIGHTED SUM.  They filter each state on
     its own criteria with AND logic; a design that fails any criterion is
     out.  A weighted composite instead lets a fatal flaw be offset by a good
     score elsewhere -- which is how 15 of the 21 round-1 G4FP variants that
     could not fold reached the order sheet.

  2. PROGRESSIVE TIGHTENING.  Their apo pLDDT threshold rises 82 -> 85 -> 86
     across rounds, with the complex states held to 92.  Cheap criteria run
     first on everything; expensive states are only predicted for survivors.

  3. THE MECHANISTICALLY ESSENTIAL SPECIES IS IN THE MODEL.  Modelling Zn(II)
     alone gave 3/89 active designs and could not distinguish active from
     inactive; adding the catalytic water gave 31/38.  For G4FP that species
     is the mature chromophore, restored by 10_graft_chromophore.py.

Their published thresholds, for reference (Supplementary Methods):

    apo protease   pLDDT > 86, catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A
    ES complex     pLDDT > 92, mean interface PAE < 5, min interface PAE < 1.2,
                   catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A
    TSA complex    pLDDT > 92, mean interface PAE < 5,
                   catalytic sidechain RMSD < 1.5 A, Ca RMSD < 1.8 A
    geometry       Zn-to-chelator distances within 0.5 A of the 4qhp crystal;
                   Zn-to-scissile-carbonyl within 0.5 A of the astacin ES model
    iteration      AF3 TSA models fed back to MPNN for one more round of design
                   and prediction; only designs passing BOTH before and after
                   were ordered

G4FP stage definitions
----------------------
The thresholds below are NOT transferred from the paper -- their numbers are
calibrated against a crystal-structure theozyme that G4FP has no equivalent
of.  They are set from the round-1 distributions and are deliberately
conservative.  Every one is a CLI flag.  Re-calibrate after round 2.

    Stage A  chromophore compatibility (cheap, structure-only, no AF3)
             mature CRO grafts without clashing; triad geometry in band
    Stage B  apo foldability
             apo pLDDT at least the parent template's (round-1 evidence:
             all 19 dark variants sat below their template, both bright
             ones above)
    Stage C  holo complex quality
             complex pLDDT, interface PAE to the G4, protein-DNA clashes
    Stage D  self-consistency
             LigandMPNN re-run on the AF3 model must recover the designed
             sequence above a threshold

Stage D requires running LigandMPNN; it is reported but not enforced unless
--sc-recovery is passed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ----------------------------------------------------------------- stages

class Stage:
    """One named gate stage: a list of (label, predicate) criteria, AND-ed."""

    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description
        self.criteria: List[tuple] = []

    def add(self, label: str, fn, required_col: Optional[str] = None) -> "Stage":
        self.criteria.append((label, fn, required_col))
        return self

    def apply(self, df: pd.DataFrame) -> tuple:
        mask = pd.Series(True, index=df.index)
        detail = []
        for label, fn, col in self.criteria:
            if col is not None and col not in df.columns:
                detail.append((label, None, "column absent -- criterion SKIPPED"))
                continue
            m = fn(df).fillna(False)
            removed = int((~m & mask).sum())
            mask = mask & m
            detail.append((label, removed, ""))
        return mask, detail


def build_cascade(a) -> List[Stage]:
    A = Stage("A_chromophore",
              "mature chromophore compatibility (structure-only, no AF3 needed)")
    A.add(f"grafted CRO min distance >= {a.cro_min_dist} A",
          lambda d: d["cro_min_dist"] >= a.cro_min_dist, "cro_min_dist")
    A.add(f"new clashes from CRO graft <= {a.cro_max_clashes}",
          lambda d: d["new_clashes"] <= a.cro_max_clashes, "new_clashes")
    A.add(f"graft fit RMSD <= {a.graft_fit_rmsd} A",
          lambda d: d["graft_fit_rmsd"] <= a.graft_fit_rmsd, "graft_fit_rmsd")

    B = Stage("B_apo_foldability",
              "the design must be able to fold and mature at all")
    B.add(f"apo pLDDT - template apo pLDDT >= {a.apo_vs_template}",
          lambda d: d["apo_plddt_vs_template"] >= a.apo_vs_template,
          "apo_plddt_vs_template")
    B.add(f"apo mean pLDDT >= {a.apo_plddt}",
          lambda d: d["apo_mean_plddt"] >= a.apo_plddt, "apo_mean_plddt")

    C = Stage("C_holo_complex", "bound-state quality and G4 interface")
    C.add(f"holo mean pLDDT >= {a.holo_plddt}",
          lambda d: d["holo_mean_plddt"] >= a.holo_plddt, "holo_mean_plddt")
    C.add(f"holo pTM >= {a.holo_ptm}", lambda d: d["holo_ptm"] >= a.holo_ptm, "holo_ptm")
    C.add(f"holo iPTM >= {a.holo_iptm}", lambda d: d["holo_iptm"] >= a.holo_iptm, "holo_iptm")
    C.add(f"mean interface PAE <= {a.interface_pae}",
          lambda d: d["holo_interface_pae"] <= a.interface_pae, "holo_interface_pae")
    C.add(f"protein-DNA clash pairs <= {a.max_clash_pairs}",
          lambda d: d["n_clash_pairs"] <= a.max_clash_pairs, "n_clash_pairs")

    D = Stage("D_self_consistency",
              "LigandMPNN re-run on the AF3 model recovers the designed sequence")
    if a.sc_recovery is not None:
        D.add(f"sequence recovery >= {a.sc_recovery}",
              lambda d: d["sc_recovery"] >= a.sc_recovery, "sc_recovery")
    return [A, B, C, D]


# ------------------------------------------------------------ stratifying

def stratify_select(df: pd.DataFrame, col: str, edges: List[float],
                    per_band: int, rank_by: str) -> pd.DataFrame:
    """Sample across bands of `col` rather than maximising it.

    Round 1 only ever tested designs BELOW their template's apo order, and
    everything that cleared it was fluorescent.  The upper edge of the viable
    window is therefore uncharacterised: too far above the floor is where
    ligand-independent folding (pKZ168) is expected.  Sampling bands locates
    that edge instead of assuming where it is.
    """
    if col not in df.columns:
        return df.nlargest(per_band * (len(edges) - 1), rank_by)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        band = df[(df[col] >= lo) & (df[col] < hi)].copy()
        band["band"] = f"[{lo},{hi})"
        out.append(band.nlargest(per_band, rank_by))
    return pd.concat(out, ignore_index=True) if out else df.head(0)


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, required=True,
                    help="rescored table from 09_rescore_and_select.py")
    ap.add_argument("--graft-report", type=Path, default=None,
                    help="10_graft_report.csv -- adds the stage-A columns")
    ap.add_argument("--out", type=Path, required=True)
    # stage A
    ap.add_argument("--cro-min-dist", type=float, default=2.3)
    ap.add_argument("--cro-max-clashes", type=int, default=0)
    ap.add_argument("--graft-fit-rmsd", type=float, default=2.0)
    # stage B
    ap.add_argument("--apo-vs-template", type=float, default=0.0)
    ap.add_argument("--apo-plddt", type=float, default=0.0)
    # stage C
    ap.add_argument("--holo-plddt", type=float, default=70.0)
    ap.add_argument("--holo-ptm", type=float, default=0.5)
    ap.add_argument("--holo-iptm", type=float, default=0.3)
    ap.add_argument("--interface-pae", type=float, default=1e9)
    ap.add_argument("--max-clash-pairs", type=int, default=10 ** 9)
    # stage D
    ap.add_argument("--sc-recovery", type=float, default=None)
    # selection
    ap.add_argument("--select", choices=["stratify", "rank"], default="stratify")
    ap.add_argument("--strat-col", default="apo_plddt_vs_template")
    ap.add_argument("--strat-edges", default="0,0.5,1,2,100")
    ap.add_argument("--per-band", type=int, default=2)
    ap.add_argument("--rank-by", default="cro_min_dist")
    ap.add_argument("--top-n", type=int, default=5, help="per template, rank mode")
    ap.add_argument("--validate", type=Path, default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.table)
    if args.graft_report and args.graft_report.exists():
        g = pd.read_csv(args.graft_report)
        g["design_id"] = g["output"].str.extract(r"_d(\d{3,5})_CRO")[0]
        g["global_id"] = g["template"] + "_design_" + g["design_id"]
        keep = ["global_id", "cro_min_dist", "new_clashes", "graft_fit_rmsd"]
        df = df.merge(g[keep], on="global_id", how="left")
        print(f"merged chromophore graft metrics for "
              f"{int(df['cro_min_dist'].notna().sum())}/{len(df)} designs")

    print(f"\n{len(df)} designs entering the cascade")
    surviving = df.copy()
    trace = []
    for st in build_cascade(args):
        if surviving.empty:
            break
        mask, detail = st.apply(surviving)
        n_before, n_after = len(surviving), int(mask.sum())
        print(f"\n[{st.name}] {st.description}")
        print(f"    {n_before} -> {n_after}")
        for label, removed, note in detail:
            print(f"      {label:52} " + (note if note else f"removed {removed}"))
        trace.append({"stage": st.name, "in": n_before, "out": n_after})
        surviving = surviving[mask].copy()

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(trace).to_csv(args.out / "11_cascade_trace.csv", index=False)
    surviving.to_csv(args.out / "11_survivors.csv", index=False)

    if surviving.empty:
        print("\nNo designs survive the cascade. Loosen the stage thresholds.")
        return

    if args.select == "stratify":
        edges = [float(x) for x in args.strat_edges.split(",")]
        sel = stratify_select(surviving, args.strat_col, edges,
                              args.per_band, args.rank_by)
        print(f"\nStratified selection on {args.strat_col}, "
              f"{args.per_band} per band, ranked within band by {args.rank_by}")
        if "band" in sel.columns:
            print(sel.groupby("band").size().to_string())
    else:
        sel = (surviving.sort_values(args.rank_by, ascending=False)
               .groupby("template").head(args.top_n))
        print(f"\nTop {args.top_n} per template by {args.rank_by}")
    sel.to_csv(args.out / "11_selected.csv", index=False)
    print(f"selected {len(sel)} designs -> {args.out}/11_selected.csv")

    if args.validate and args.validate.exists():
        v = pd.read_csv(args.validate)
        norm = lambda g: str(g)[5:] if str(g).startswith("G4FP_") else str(g)
        v["_k"] = v["global_id"].map(norm); df["_k"] = df["global_id"].map(norm)
        m = df.merge(v.drop(columns="global_id"), on="_k", how="inner")
        if not m.empty and args.rank_by in m.columns:
            m = m.sort_values(args.rank_by, ascending=False).reset_index(drop=True)
            m["rank"] = m.index + 1
            cols = [c for c in ("rank", "global_id", "outcome", args.rank_by,
                                "cro_min_dist", "apo_plddt_vs_template")
                    if c in m.columns]
            print(f"\nValidation on {len(m)} tested variants, ranked by {args.rank_by}:")
            print(m[cols].to_string(index=False))
            m.to_csv(args.out / "11_validation.csv", index=False)


if __name__ == "__main__":
    main()

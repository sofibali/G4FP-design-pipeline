#!/usr/bin/env python3
"""
12_iterate.py -- iterative filtering pipeline for G4FP design.

Goal being optimised
--------------------
A reporter that is fluorescent ONLY when bound to a G-quadruplex.  Three
properties have to hold simultaneously, and they are tested as separate gate
groups rather than traded off against each other in one score:

  FOLD   the holo state is a well-formed, maturable beta-barrel
  BIND   the protein actually engages the G4
  SWITCH the chromophore environment differs substantially between holo and apo
  APO    the apo state sits in a window: ordered enough that the chromophore
         can mature at all, not so ordered that fluorescence is
         ligand-independent

Why the switch gate is chromophore-local, not global
----------------------------------------------------
Measured across the 1,000 round-1 designs, global unfolding fights
maturation but local chromophore rearrangement does not:

    metric                    corr with apo_plddt_vs_template   clear apo floor
                                                                 AND top quartile
    mean_plddt_diff                   -0.321 (p=2e-25)                20/280
    global_rmsd                       -0.052 (p=0.10)                 44/280
    chromophore_plddt_diff            +0.277 (p=4e-19)               106/280
    chromophore_rmsd                  +0.518 (p=1e-69)               145/280

A design cannot both unfold globally in the apo state and still mature its
chromophore -- maturation requires the folded barrel.  Requiring a large
*global* holo-apo difference therefore selects against the very thing that
makes the protein fluorescent, and leaves almost nothing (20 of 1,000).  The
chromophore-local measures carry the switch without that conflict, which is
also the GCaMP mechanism: the barrel stays, the chromophore environment
changes.

Iteration
---------
Modelled on Chen, Wu et al. (Nature 2026) which iterates sequence design
against structure prediction with thresholds tightening each round --
apo pLDDT > 82, then > 85, then > 86 -- and requires the final designs to pass
BOTH BEFORE AND AFTER a last design/prediction cycle.

One G4FP iteration is:

    survivors -> LigandMPNN redesign on the AF3 model (10_graft first)
              -> AF3 holo + apo (lab server)
              -> 09_rescore_and_select.py
              -> this script, one notch tighter
              -> survivors

Thresholds are PERCENTILES of the current pool, not absolute numbers, with a
small set of absolute floors for physical plausibility.  The paper's absolute
cutoffs are calibrated against a crystal-structure theozyme and do not
transfer; percentiles self-calibrate to whatever the corrected metrics look
like, and tighten monotonically by construction.

Usage
-----
    # first pass over a freshly rescored pool
    python 12_iterate.py --table 09_results_full/09_all_designs_rescored.csv \\
                         --iteration 1 --out iter1

    # later rounds: --previous enforces "passed before AND after"
    python 12_iterate.py --table iter2_rescored.csv --iteration 2 \\
                         --previous iter1/12_survivors.csv --out iter2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Threshold ladder
# --------------------------------------------------------------------------
# Each entry: fraction of the pool to KEEP on that criterion at that iteration.
# Tightens monotonically.  Absolute floors are applied in addition and never
# relax.  Edit here, or override any single value on the command line.

LADDER: Dict[int, Dict[str, float]] = {
    1: {"fold": 0.60, "bind": 0.60, "switch": 0.50, "apo_hi": 0.95},
    2: {"fold": 0.40, "bind": 0.40, "switch": 0.35, "apo_hi": 0.90},
    3: {"fold": 0.25, "bind": 0.25, "switch": 0.25, "apo_hi": 0.85},
}

ABSOLUTE = {
    "holo_mean_plddt_min": 70.0,   # must be a plausible fold when bound
    "holo_ptm_min": 0.50,
    "holo_iptm_min": 0.30,
    "apo_vs_template_min": 0.0,    # lower edge of the apo window: below this,
                                   # all 19 experimentally dark round-1
                                   # variants; see README section 2
}

# Column fallbacks: corrected names first, round-1 names second.
COLS = {
    "holo_plddt":  ["holo_mean_plddt"],
    "holo_ptm":    ["holo_ptm"],
    "holo_iptm":   ["holo_iptm"],
    "holo_chrom":  ["holo_chromophore_plddt_fixed", "holo_chromophore_plddt"],
    "iface_pae":   ["holo_interface_pae"],
    "apo_rel":     ["apo_plddt_vs_template"],
    "chrom_diff":  ["chromophore_plddt_diff_fixed", "chromophore_plddt_diff"],
    "chrom_rmsd":  ["chromophore_rmsd"],
}


def pick(df: pd.DataFrame, key: str) -> Optional[str]:
    for c in COLS[key]:
        if c in df.columns and df[c].notna().any():
            return c
    return None


def pct_gate(df: pd.DataFrame, col: str, keep: float, higher_is_better: bool):
    """Boolean mask keeping the best `keep` fraction on `col`."""
    v = df[col]
    if v.notna().sum() == 0:
        return None, None
    thr = v.quantile(1 - keep) if higher_is_better else v.quantile(keep)
    return (v >= thr) if higher_is_better else (v <= thr), float(thr)


# --------------------------------------------------------------------------
# Gate groups
# --------------------------------------------------------------------------

def run_groups(df: pd.DataFrame, keep: Dict[str, float], absolute: Dict[str, float],
               apo_hi_abs: Optional[float]) -> (pd.DataFrame, List[dict]):
    trace: List[dict] = []
    mask = pd.Series(True, index=df.index)

    def apply(group: str, label: str, m, thr=None):
        nonlocal mask
        if m is None:
            trace.append({"group": group, "criterion": label, "threshold": None,
                          "removed": None, "note": "metric unavailable -- SKIPPED"})
            return
        removed = int((~m & mask).sum())
        mask = mask & m.fillna(False)
        trace.append({"group": group, "criterion": label,
                      "threshold": None if thr is None else round(thr, 3),
                      "removed": removed, "n_remaining": int(mask.sum()), "note": ""})

    # ---- APO: a window, not a floor ------------------------------------
    c = pick(df, "apo_rel")
    if c:
        apply("APO", f"{c} >= {absolute['apo_vs_template_min']} (lower edge: "
                     "below this sat all 19 dark round-1 variants)",
              df[c] >= absolute["apo_vs_template_min"], absolute["apo_vs_template_min"])
        if apo_hi_abs is not None:
            apply("APO", f"{c} <= {apo_hi_abs} (upper edge, explicit)",
                  df[c] <= apo_hi_abs, apo_hi_abs)
        else:
            m, thr = pct_gate(df, c, keep["apo_hi"], higher_is_better=False)
            apply("APO", f"{c} <= P{int(keep['apo_hi']*100)} (upper edge: too far "
                         "above template is ligand-independent folding)", m, thr)

    # ---- FOLD: the ON state must be a real barrel ----------------------
    for key, lab, absmin in (("holo_plddt", "holo mean pLDDT", absolute["holo_mean_plddt_min"]),
                             ("holo_ptm", "holo pTM", absolute["holo_ptm_min"])):
        c = pick(df, key)
        if c:
            apply("FOLD", f"{c} >= {absmin} (absolute)", df[c] >= absmin, absmin)
            m, thr = pct_gate(df, c, keep["fold"], True)
            apply("FOLD", f"{c} >= P{int((1-keep['fold'])*100)}", m, thr)
    c = pick(df, "holo_chrom")
    if c:
        m, thr = pct_gate(df, c, keep["fold"], True)
        apply("FOLD", f"{c} >= P{int((1-keep['fold'])*100)} (chromophore ordered "
                      "in the ON state)", m, thr)

    # ---- BIND: the protein must engage the G4 --------------------------
    c = pick(df, "holo_iptm")
    if c:
        apply("BIND", f"{c} >= {absolute['holo_iptm_min']} (absolute)",
              df[c] >= absolute["holo_iptm_min"], absolute["holo_iptm_min"])
        m, thr = pct_gate(df, c, keep["bind"], True)
        apply("BIND", f"{c} >= P{int((1-keep['bind'])*100)}", m, thr)
    c = pick(df, "iface_pae")
    if c:
        m, thr = pct_gate(df, c, keep["bind"], False)
        apply("BIND", f"{c} <= P{int(keep['bind']*100)} (interface error)", m, thr)

    # ---- SWITCH: chromophore environment must change on binding --------
    for key, lab in (("chrom_diff", "chromophore pLDDT holo-apo"),
                     ("chrom_rmsd", "chromophore RMSD holo vs apo")):
        c = pick(df, key)
        if c:
            m, thr = pct_gate(df, c, keep["switch"], True)
            apply("SWITCH", f"{c} >= P{int((1-keep['switch'])*100)} ({lab})", m, thr)

    return df[mask].copy(), trace


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, required=True,
                    help="rescored metric table (output of 09_rescore_and_select.py)")
    ap.add_argument("--iteration", type=int, default=1, choices=sorted(LADDER))
    ap.add_argument("--previous", type=Path, default=None,
                    help="survivors of the previous iteration; enforces the "
                         "'passed before AND after' rule")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep-fold", type=float, default=None)
    ap.add_argument("--keep-bind", type=float, default=None)
    ap.add_argument("--keep-switch", type=float, default=None)
    ap.add_argument("--apo-hi", type=float, default=None,
                    help="explicit upper edge for apo_plddt_vs_template; "
                         "default is a percentile of the pool")
    ap.add_argument("--validate", type=Path, default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.table)
    keep = dict(LADDER[args.iteration])
    for k in ("fold", "bind", "switch"):
        v = getattr(args, f"keep_{k}")
        if v is not None:
            keep[k] = v
    print(f"iteration {args.iteration}: {len(df)} designs in, keep fractions {keep}")

    missing = [k for k in COLS if pick(df, k) is None]
    if missing:
        print(f"  metrics unavailable (their criteria are skipped): {', '.join(missing)}")

    surv, trace = run_groups(df, keep, ABSOLUTE, args.apo_hi)

    last = None
    for t in trace:
        if t["group"] != last:
            print(f"\n[{t['group']}]")
            last = t["group"]
        if t["note"]:
            print(f"    {t['criterion'][:70]:72} {t['note']}")
        else:
            print(f"    {t['criterion'][:70]:72} thr={t['threshold']} "
                  f"removed {t['removed']} -> {t['n_remaining']}")
    print(f"\nsurvivors: {len(surv)} of {len(df)}")

    if args.previous and args.previous.exists():
        prev = pd.read_csv(args.previous)
        before = set(prev["global_id"])
        surv["passed_previous_iteration"] = surv["global_id"].isin(before)
        n_both = int(surv["passed_previous_iteration"].sum())
        print(f"passed BOTH this and the previous iteration: {n_both}")
        surv = surv[surv["passed_previous_iteration"]].copy()

    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(trace).to_csv(args.out / "12_gate_trace.csv", index=False)
    surv.to_csv(args.out / "12_survivors.csv", index=False)
    surv[["global_id"]].to_csv(args.out / "12_next_round_inputs.csv", index=False)
    json.dump({"iteration": args.iteration, "keep": keep, "absolute": ABSOLUTE,
               "n_in": len(df), "n_out": len(surv)},
              open(args.out / "12_config.json", "w"), indent=2)

    if "template" in surv.columns and len(surv):
        print("\nsurvivors per template:")
        print(surv.groupby("template").size().to_string())

    if args.validate and args.validate.exists():
        v = pd.read_csv(args.validate)
        norm = lambda g: str(g)[5:] if str(g).startswith("G4FP_") else str(g)
        v["_k"] = v["global_id"].map(norm)
        allk = df.assign(_k=df["global_id"].map(norm)).merge(
            v.drop(columns="global_id"), on="_k", how="inner")
        survk = set(surv["global_id"].map(norm))
        allk["survived"] = allk["_k"].isin(survk)
        print(f"\nvalidation on {len(allk)} experimentally tested variants:")
        print(allk.groupby(["outcome", "survived"]).size().to_string())
        allk.to_csv(args.out / "12_validation.csv", index=False)

    print(f"\nwrote {args.out}/ -- feed 12_next_round_inputs.csv to "
          f"10_graft_chromophore.py --only for the next design round")


if __name__ == "__main__":
    main()

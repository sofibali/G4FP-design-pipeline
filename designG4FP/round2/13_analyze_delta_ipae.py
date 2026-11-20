#!/usr/bin/env python3
"""
13_analyze_delta_ipae.py -- what the confidence-contrast gate says about round 1.

Reads 13_delta_ipae.csv (from 13_delta_ipae_gate.py) and tested_variants.csv,
and answers three questions in order.  The order matters: if the first answer
is no, the other two are not worth acting on.

  Q1  Does the round-1 population have ANY dPAE spread?
      If every design scores ~0, the architecture does not switch and
      re-gating it cannot find a switch that is not there.

  Q2  Where do the 21 experimentally tested designs sit in that spread?

  Q3  Could this gate have been validated on round-1 data at all?
      It could not.  The gate predicts SWITCHING.  Of the 21 tested, 2 were
      fluorescent and 0 were G4-responsive, so the switch phenotype has zero
      positives.  The script states the bound rather than computing a p-value
      against a label that does not exist.

Usage:  python3 13_analyze_delta_ipae.py
"""

from __future__ import annotations

import csv
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GATE_CSV = os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1 else "13_delta_ipae_module.csv")
TESTED_CSV = os.path.join(HERE, "tested_variants.csv")


def read_gate(path):
    rows = []
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    for r in csv.DictReader(lines):
        # output_G4FP_x + design_7  ->  G4FP_x_design_0007
        tpl = r["template"][len("output_"):] if r["template"].startswith("output_") else r["template"]
        n = int(r["design_id"].rsplit("_", 1)[1])
        r["global_id"] = "%s_design_%04d" % (tpl, n)
        for k in ("ipae_apo", "ipae_holo", "dpae", "dpae_contact", "bind_ipae_holo",
                  "bind_pae_min", "iptm_holo", "n_module_contact"):
            r[k] = float(r[k]) if r[k] not in ("", None) else float("nan")
        rows.append(r)
    return rows


def read_tested(path):
    with open(path) as fh:
        return {r["global_id"]: r["outcome"] for r in csv.DictReader(fh)}


def pct(sorted_vals, v):
    """Percentile of v within sorted_vals, 0-100."""
    lo = 0
    for x in sorted_vals:
        if x < v:
            lo += 1
        else:
            break
    return 100.0 * lo / len(sorted_vals)


def describe(name, vals, unit="A"):
    vals = [v for v in vals if v == v]
    if not vals:
        print("  %-16s no values" % name)
        return
    s = sorted(vals)
    q = lambda p: s[min(len(s) - 1, int(p * len(s)))]
    print("  %-16s n=%4d  min %7.3f  p5 %7.3f  median %7.3f  p95 %7.3f  max %7.3f  sd %6.3f  (%s)"
          % (name, len(s), s[0], q(.05), q(.50), q(.95), s[-1], statistics.pstdev(s), unit))


def main():
    if not os.path.exists(GATE_CSV):
        print("missing %s -- run 13_delta_ipae_gate.py first" % GATE_CSV, file=sys.stderr)
        return 1
    rows = read_gate(GATE_CSV)
    tested = read_tested(TESTED_CSV)

    print("=" * 78)
    print("Q1  Is there any dPAE spread in the round-1 population?")
    print("=" * 78)
    describe("dpae (window)", [r["dpae"] for r in rows])
    describe("dpae (contact)", [r["dpae_contact"] for r in rows])
    describe("ipae_apo", [r["ipae_apo"] for r in rows])
    describe("ipae_holo", [r["ipae_holo"] for r in rows])
    describe("bind_ipae_holo", [r["bind_ipae_holo"] for r in rows])
    describe("n_module_contact", [r["n_module_contact"] for r in rows], unit="residues")

    dp = sorted(r["dpae"] for r in rows if r["dpae"] == r["dpae"])
    n_hits = sum(1 for v in dp if v > 3.0)
    print("\n  Chang & Polizzi's published switch threshold is dPAE > 3.")
    print("  Designs clearing it here: %d / %d  (%.1f%%)" % (n_hits, len(dp), 100.0 * n_hits / len(dp)))
    print("  Interquartile width of dpae: %.3f A" % (dp[int(.75 * len(dp))] - dp[int(.25 * len(dp))]))

    print()
    print("=" * 78)
    print("Q2  Where do the 21 tested designs sit?")
    print("=" * 78)
    by_id = {r["global_id"]: r for r in rows}
    found = [(gid, out, by_id[gid]) for gid, out in tested.items() if gid in by_id]
    missing = [gid for gid in tested if gid not in by_id]
    print("  matched %d of %d tested designs%s"
          % (len(found), len(tested), ("  (missing: %s)" % ", ".join(missing)) if missing else ""))

    s_dp = sorted(r["dpae"] for r in rows if r["dpae"] == r["dpae"])
    s_bind = sorted(r["bind_ipae_holo"] for r in rows if r["bind_ipae_holo"] == r["bind_ipae_holo"])
    print("\n  %-34s %-18s %8s %7s %9s %7s" % ("global_id", "outcome", "dpae", "pctile", "bind_ipae", "pctile"))
    for gid, out, r in sorted(found, key=lambda t: -t[2]["dpae"]):
        print("  %-34s %-18s %8.3f %6.1f%% %9.3f %6.1f%%"
              % (gid, out, r["dpae"], pct(s_dp, r["dpae"]),
                 r["bind_ipae_holo"], pct(s_bind, r["bind_ipae_holo"])))

    bright = [r for gid, out, r in found if out.startswith("fluorescent")]
    dark = [r for gid, out, r in found if out == "dark"]
    if bright and dark:
        print("\n  fluorescent mean dpae %.3f   dark mean dpae %.3f   difference %+.3f"
              % (statistics.mean(x["dpae"] for x in bright),
                 statistics.mean(x["dpae"] for x in dark),
                 statistics.mean(x["dpae"] for x in bright) - statistics.mean(x["dpae"] for x in dark)))

    print()
    print("=" * 78)
    print("Q3  Can this gate be validated on round-1 data?")
    print("=" * 78)
    n, k = len(found), len(bright)
    floor = 1.0 / math.comb(n, k) if 0 < k < n else float("nan")
    print("  The gate predicts G4-driven SWITCHING.")
    print("  Round-1 outcomes: %d tested, %d fluorescent, 0 G4-responsive." % (n, k))
    print("  There are no positives for the phenotype this gate predicts, so no")
    print("  validation against round-1 data is possible -- not a weak one, none.")
    print()
    print("  Even for the FLUORESCENCE label, an exact permutation test over")
    print("  C(%d,%d) = %d arrangements has a floor of p = %.4f." % (n, k, math.comb(n, k), floor))
    print("  No metric can reach significance at this sample size; a p-value here")
    print("  bounds a claim, it does not support one.")
    print()
    print("  => Adopt or reject the gate on mechanism and on Chang & Polizzi's")
    print("     prospective result, then validate it in round 2 by testing designs")
    print("     chosen BECAUSE of it. Scoring the round-1 population answers only")
    print("     Q1: whether this architecture has a switch to find.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

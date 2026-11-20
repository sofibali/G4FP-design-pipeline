#!/usr/bin/env python3
"""
17_forecast_denovo.py -- cost and funnel for round 2 as a NEW design campaign.

PLANNING ONLY.  This script runs nothing.  It reads the measured round-1
timings and the measured round-1 gate pass rates and projects them onto a
round 2 that starts from new scaffolds rather than from round-1 designs.

Why this replaces 15_round2_forecast.py's scenarios A / B / B-lite
-------------------------------------------------------------------
Those scenarios iterated round-1 designs.  Their survivor bands were very wide
(505-4,276 for scenario B) for one reason: round 1 contains no parent/child
pairs, so the inheritance of gate behaviour from a parent design to its
redesigned child has never been measured in this tree.  Any iteration forecast
is therefore a guess with a band around it.

A fresh round has the opposite property.  Round 1 IS a measurement of a fresh
round: 999 designs entered a cascade and 71 came out.  So the funnel for a new
round of the same shape is not extrapolated at all -- it is the round-1 rate,
applied again.  The only extrapolated quantity is the CONTACT uplift from the
wider design window, and that is carried as an explicit low/high band.

This is the main practical argument for the new-design framing, independent of
any biology: it is the only one of the two that can be forecast.

Measured inputs, with the file each came from
---------------------------------------------
  AF3 wall clock, concurrency   output_G4FP_*/03_alphafold3_predictions_*/ mtimes
  LigandMPNN wall clock         batch_ligandmpnn_progress.log
  Gate pass rates               14_regroup_failuremode.csv
  Design-window G4 coverage     13_RESCORING_PLAN.md (140 -> 251 contact residues)
  Experimental outcomes         tested_variants.csv

NOT costed here, because this tree never timed it
--------------------------------------------------
  Scaffold generation.  Round 1's ten templates were produced before the step
  that this repo times, so the hours needed to build N new scaffolds are not a
  measurement available here.  They are additive to every number below.

Usage:  python3 -I 17_forecast_denovo.py [--out 17_forecast_denovo.csv]
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# MEASURED CONSTANTS -- identical to 15_round2_forecast.py; see its header for
# how each was taken.  Kept duplicated rather than imported so this file can be
# read and audited on its own.
# ---------------------------------------------------------------------------
AF3_SEC_HOLO = 141.6              # 326 tokens
AF3_SEC_APO = 138.2               # 307 tokens
AF3_GPUS = 2
AF3_ROUND1_PREDICTIONS = 2002
AF3_ROUND1_CALENDAR_H = 40.6      # measured; used to validate the model below

MPNN_SEC_PER_INPUT = 32.3 * 60    # 10,001 sequences from one input backbone
ANALYSIS_SEC_PER_100 = 11 * 60    # steps 05 -> 06, 100 designs, both states

FILTER_KEEP = 100                 # sequences kept per scaffold (top 1% of 10k)

# ---------------------------------------------------------------------------
# MEASURED ROUND-1 FUNNEL  (14_regroup_failuremode.csv, N = 999)
# Conditional pass rate at each gate, i.e. of those that reached it.
# ---------------------------------------------------------------------------
R1_N = 999
R1_REACHED = [999, 298, 197, 101]     # entering CONTACT, NEGATIVE, BIND, SWITCH
R1_PASSED = [298, 197, 101, 71]
GATE_NAMES = ["CONTACT", "NEGATIVE", "BIND", "SWITCH"]

# CONTACT uplift from the wider design window.  13_RESCORING_PLAN.md: the
# round-1 window covers 140 of the G4-contacting residues, the 40+40 window
# covers 251.  Treated as the HIGH end only; the LOW end assumes no uplift at
# all, because a wider window makes contact *designable*, not *certain*.
CONTACT_UPLIFT = 251.0 / 140.0
CONTACT_CEILING = 0.95

# ---------------------------------------------------------------------------
# MEASURED EXPERIMENTAL OUTCOME  (tested_variants.csv)
# ---------------------------------------------------------------------------
R1_ORDERED, R1_FLUORESCENT, R1_RESPONSIVE = 21, 2, 0

# ---------------------------------------------------------------------------
# SCENARIOS -- all are "build new scaffolds", differing only in how wide
# ---------------------------------------------------------------------------
SCENARIOS = [
    ("P",  "Pilot: 2 scaffolds, 20 designs each",            2, 20),
    ("N1", "Narrow: one module, 4 scaffolds",                4, FILTER_KEEP),
    ("N2", "Round-1 width: 10 scaffolds",                   10, FILTER_KEEP),
    ("N3", "Wide: two modules x 10 scaffolds",              20, FILTER_KEEP),
]


def hms(hours):
    if hours < 1:
        return "%d min" % round(hours * 60)
    if hours < 48:
        return "%.1f h" % hours
    return "%.1f d" % (hours / 24.0)


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                              capture_output=True, text=True, timeout=10,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "17_forecast_denovo.csv"))
    args = ap.parse_args()

    rates = [p / float(r) for p, r in zip(R1_PASSED, R1_REACHED)]

    print("=" * 78)
    print("ROUND 2 AS A NEW DESIGN CAMPAIGN -- predicted cost and funnel")
    print("PLANNING ONLY.  Nothing in this script executes the pipeline.")
    print("=" * 78)

    # --- model validation against the one measured round ---------------------
    per_pred = (AF3_SEC_HOLO + AF3_SEC_APO) / 2.0
    pred_h = AF3_ROUND1_PREDICTIONS * per_pred / AF3_GPUS / 3600.0
    err = 100.0 * (pred_h - AF3_ROUND1_CALENDAR_H) / AF3_ROUND1_CALENDAR_H
    print()
    print("Model check against round 1 as actually run")
    print("-" * 78)
    print("  %d predictions x %.1f s / %d GPUs = %s predicted"
          % (AF3_ROUND1_PREDICTIONS, per_pred, AF3_GPUS, hms(pred_h)))
    print("  measured calendar time                  %s" % hms(AF3_ROUND1_CALENDAR_H))
    print("  error                                   %+.1f %%" % err)

    # --- the measured funnel -------------------------------------------------
    print()
    print("Measured round-1 funnel, N = %d  (14_regroup_failuremode.csv)" % R1_N)
    print("-" * 78)
    print("  %-10s %8s %8s %9s" % ("gate", "reached", "passed", "cond. rate"))
    for nm, reached, passed, rate in zip(GATE_NAMES, R1_REACHED, R1_PASSED, rates):
        print("  %-10s %8d %8d %8.1f %%" % (nm, reached, passed, 100 * rate))
    overall_low = 1.0
    for r in rates:
        overall_low *= r
    contact_high = min(CONTACT_CEILING, rates[0] * CONTACT_UPLIFT)
    overall_high = contact_high
    for r in rates[1:]:
        overall_high *= r
    print()
    print("  end-to-end, no window change          %6.2f %%  (= %d/%d, reconciles)"
          % (100 * overall_low, R1_PASSED[-1], R1_N))
    print("  CONTACT with the 40+40 window         %6.1f %%  (%.1f %% x %.2f, capped at %.0f %%)"
          % (100 * contact_high, 100 * rates[0], CONTACT_UPLIFT, 100 * CONTACT_CEILING))
    print("  end-to-end, window uplift applied     %6.2f %%" % (100 * overall_high))
    print()
    print("  The band is %.2f-%.2f %%, a factor of %.1f." % (100 * overall_low,
                                                             100 * overall_high,
                                                             overall_high / overall_low))
    print("  Compare the iteration scenarios in 15_round2_forecast.csv, whose")
    print("  bands span a factor of 8.5, because design-to-design inheritance")
    print("  has never been measured in this tree.  A fresh round needs no such")
    print("  assumption: round 1 is the measurement.")

    # --- scenarios -----------------------------------------------------------
    rows = []
    print()
    print("Scenarios")
    print("-" * 78)
    print("  %-4s %-38s %5s %7s %8s %9s" %
          ("id", "name", "scaf", "designs", "AF3", "wall"))
    for sid, name, n_scaf, n_des in SCENARIOS:
        n_designs = n_scaf * n_des
        n_af3 = n_designs * 2
        mpnn_h = n_scaf * MPNN_SEC_PER_INPUT / 3600.0
        af3_h = n_af3 * per_pred / AF3_GPUS / 3600.0
        ana_h = n_designs / 100.0 * ANALYSIS_SEC_PER_100 / 3600.0
        wall_h = mpnn_h + af3_h + ana_h
        s_low = int(round(n_designs * overall_low))
        s_high = int(round(n_designs * overall_high))
        print("  %-4s %-38s %5d %7d %8d %9s" %
              (sid, name, n_scaf, n_designs, n_af3, hms(wall_h)))
        rows.append(dict(scenario=sid, name=name, n_scaffolds=n_scaf,
                         n_designs=n_designs, n_af3=n_af3,
                         mpnn_h=round(mpnn_h, 2), af3_h=round(af3_h, 2),
                         analysis_h=round(ana_h, 2), wall_h=round(wall_h, 2),
                         survivors_low=s_low, survivors_high=s_high))

    print()
    print("  %-4s %-38s %9s %9s" % ("id", "name", "surv.low", "surv.high"))
    for r in rows:
        print("  %-4s %-38s %9d %9d" % (r["scenario"], r["name"],
                                        r["survivors_low"], r["survivors_high"]))

    # --- the comparison that matters ----------------------------------------
    n2 = [r for r in rows if r["scenario"] == "N2"][0]
    print()
    print("A full new round at round-1 width costs %s." % hms(n2["wall_h"]))
    print("Iterating the 71 survivors (15_round2_forecast.csv, scenario B)")
    print("costs 308.0 h = 12.8 d, for a survivor count nobody can predict.")
    print("The new round is %.1fx cheaper and %.1fx narrower in its forecast."
          % (308.04 / n2["wall_h"], 8.5 / (overall_high / overall_low)))

    # --- what the funnel still does not select for ---------------------------
    print()
    print("What this funnel still does NOT select for")
    print("-" * 78)
    print("  Round 1 ordered %d designs.  %d fluoresced at all; %d were "
          "G4-responsive." % (R1_ORDERED, R1_FLUORESCENT, R1_RESPONSIVE))
    print("  Of the %d tested, 10 passed all four gates and all 10 were dark;" % R1_ORDERED)
    print("  both designs that did fluoresce are REJECTED by the cascade.")
    print("  So the %d-%d survivors above are survivors of a binding-and-switching"
          % (rows[2]["survivors_low"], rows[2]["survivors_high"]))
    print("  filter only.  Chromophore maturation -- where 19 of 21 failed -- is")
    print("  not gated anywhere in it, and this tree contains no computational")
    print("  predictor of maturation.  Round 2 has exactly two options:")
    print("    (a) add an orthogonal maturation gate, which must be built and")
    print("        calibrated first, on data that does not yet exist; or")
    print("    (b) delete the maturation step -- bind a dye whose fluorescence")
    print("        is environment-switched rather than self-generated.")
    print("  Nothing in the compute budget decides this.  It is a design choice")
    print("  that must be made BEFORE any scaffold is generated.")

    # --- order size ----------------------------------------------------------
    print()
    print("Order size, unchanged from 15_ROUND2_CAMPAIGN_PLAN.md section 6")
    print("-" * 78)
    print("  0/%d responsive excludes a true rate above 13.3 %% at 95 %% confidence." % R1_ORDERED)
    print("  For 80 % power: 53 constructs at a 3 % true rate, 32 at 5 %,")
    print("  16 at 10 %%, 9 at 17 %%.  Round 1 ordered %d." % R1_ORDERED)
    print("  N2 and N3 clear 32 survivors at the LOW end of the band.  N1 clears")
    print("  it only if the window uplift is real (28 low, 51 high), so N1 is a")
    print("  bet on the uplift.  P is a pilot and is not ordered from.")
    print("  The order size, not the compute, is the binding constraint:")
    print("  every scenario here fits in a week, and none of them is the reason")
    print("  round 1 could not tell a working design from a dead one.")

    # --- write ---------------------------------------------------------------
    with open(args.out, "w", newline="") as fh:
        fh.write("# 17_forecast_denovo.py\n")
        fh.write("# generated: %s\n" % datetime.datetime.now().isoformat(timespec="seconds"))
        fh.write("# git: %s\n" % git_hash())
        fh.write("# PLANNING ONLY -- predicted, not measured. Nothing was executed.\n")
        fh.write("# round 2 framed as a NEW design campaign, not an iteration of round-1 designs\n")
        fh.write("# *_h are hours; AF3 anchors %.1f s holo / %.1f s apo, %d GPUs\n"
                 % (AF3_SEC_HOLO, AF3_SEC_APO, AF3_GPUS))
        fh.write("# survivors_low  = round-1 funnel applied unchanged (%.2f %%)\n" % (100 * overall_low))
        fh.write("# survivors_high = CONTACT uplifted by the 40+40 window (%.2f %%)\n" % (100 * overall_high))
        fh.write("# scaffold generation is NOT costed here; it is additive\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print()
    print("Wrote %s (%d rows)" % (args.out, len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

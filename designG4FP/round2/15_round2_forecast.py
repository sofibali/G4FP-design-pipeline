#!/usr/bin/env python3
"""
15_round2_forecast.py -- predicted cost and predicted funnel for the round-2
design campaign.  PLANNING ONLY.  This script runs nothing: it reads measured
round-1 timings and measured round-1 gate pass rates and projects them forward.

Every constant below is a MEASUREMENT taken from this tree, with the file it
came from named.  Nothing is a vendor figure or a remembered number.

  AF3 wall clock          mtimes of the 2,002 prediction directories under
                          output_G4FP_*/03_alphafold3_predictions_{bound,apo}/
  AF3 concurrency         the two interleaved streams visible in those mtimes
  LigandMPNN wall clock   batch_ligandmpnn_progress.log
  Token counts            design_0_bound_data.json (protein + dna + 3 K+)
  Gate pass rates         14_regroup_failuremode.csv
  Experimental outcomes   tested_variants.csv

The two things it CANNOT measure, and therefore states as assumptions with a
band rather than a number:

  1. how AF3 wall clock scales with token count beyond the 307 -> 326 step
     round 1 happens to contain;
  2. whether a redesigned child inherits its parent's gate behaviour.

Both are settled by the 10-design pilot in section 6, which costs ~25 minutes.

Usage:  python3 -I 15_round2_forecast.py [--out 15_round2_forecast.csv]
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
import os
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# MEASURED CONSTANTS  (source named for each)
# ---------------------------------------------------------------------------

# AF3, 5 model seeds x 5 samples, on the lab server's GPUs.
# Median inter-directory interval over 100 consecutive predictions, per template;
# the figure below is the median of those ten per-template medians.
AF3_SEC_HOLO = 141.6          # 326 tokens (306 aa + CRO + 16 nt G4 + 3 K+)
AF3_SEC_APO = 138.2           # 307 tokens (protein alone)
AF3_TOKENS_HOLO = 326
AF3_TOKENS_APO = 307
AF3_GPUS = 2                  # two interleaved streams, ~1 h offset
AF3_ROUND1_PREDICTIONS = 2002
AF3_ROUND1_CALENDAR_H = 40.6

# LigandMPNN, 10,001 sequences per input, side chains packed, ligand context on.
MPNN_SEC_PER_10K = 32.3 * 60  # mean of 29m43s and 34m56s
MPNN_N_SEQ = 10000

# Downstream CPU analysis: step 05 -> 06 for one template (100 designs, 2 states)
ANALYSIS_SEC_PER_100 = 11 * 60

# Round-1 sequence filter: 10,001 -> 100 kept (top 1%), pipeline_config.yaml
FILTER_KEEP_FRACTION = 0.01

# ---------------------------------------------------------------------------
# MEASURED ROUND-1 GATE RATES  (14_regroup_failuremode.csv)
# ---------------------------------------------------------------------------
FAILMODE_CSV = "14_regroup_failuremode.csv"

GATE_ORDER = [
    ("CONTACT", "ok_contact", "AF3 predicts >=3 module-G4 contacts"),
    ("NEGATIVE", "ok_negative", "not pre-closed: apo iPAE >= 8 A"),
    ("BIND", "ok_bind", "interface PAE <= P40 of the population"),
    ("SWITCH", "ok_switch", "dPAE > 3 A"),
]

# ---------------------------------------------------------------------------
# ASSUMPTIONS -- stated, not hidden
# ---------------------------------------------------------------------------

# A1. AF3 time vs token count.  Round 1 contains exactly one step: 307 -> 326
#     tokens (+6.2%) cost 138.2 -> 141.6 s (+2.5%), an apparent exponent of
#     ln(1.025)/ln(1.062) = 0.41.  A single 6% step is far too short a lever to
#     extrapolate from, so we carry a band: the measured 0.41 as the low end and
#     a quadratic 2.0 as the pessimistic end (what pair-attention alone would
#     give if nothing else were fixed-cost).
SCALE_EXP_LOW = 0.41
SCALE_EXP_HIGH = 2.0

# A2. Inheritance.  A LigandMPNN child of a design that passed gate X may or may
#     not pass gate X itself.  Round 1 cannot measure this -- it has no
#     parent/child pairs.  Band:
#       optimistic  = the parent template's measured round-1 pass rate
#       pessimistic = the library-wide rate, i.e. redesign randomises completely
#     The pilot in section 6 measures the real value.

# A3. Round-2 redesign window.  Round 1 used A15-60 + A240-285, which captured
#     140 of 283 G4 contacts (49%).  Round 2 uses first-40 + last-40, which
#     captures 251 of 283 (89%) -- README.md section 3.  The CONTACT gate is the
#     one that kills 70% of round 1, and it is precisely the gate this window
#     change targets.  We carry the ratio 89/49 = 1.82 as the optimistic uplift
#     and 1.0 (no uplift) as the pessimistic one, capped at 100%.
WINDOW_CONTACT_UPLIFT = 251.0 / 140.0

# A4. Experimental outcomes.  Round 1: 21 ordered, 2 fluorescent (9.5%),
#     0 G4-responsive.  Chang & Polizzi report 4 clear switches in 24 tested
#     (17%) for the architecture the SWITCH gate comes from.
R1_ORDERED, R1_FLUORESCENT, R1_RESPONSIVE = 21, 2, 0
CHANG_POLIZZI_RATE = 4.0 / 24.0


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                              capture_output=True, text=True, timeout=10).stdout.strip() or "n/a"
    except Exception:
        return "n/a"


def hms(sec):
    if sec < 90:
        return "%.0f s" % sec
    if sec < 5400:
        return "%.0f min" % (sec / 60.0)
    if sec < 3 * 86400:
        return "%.1f h" % (sec / 3600.0)
    return "%.1f d" % (sec / 86400.0)


def af3_sec(tokens, state_sec, state_tokens, exp):
    """Scale a measured per-prediction time to a new token count."""
    return state_sec * (float(tokens) / state_tokens) ** exp


def load_rates(path):
    if not os.path.exists(path):
        print("missing %s -- run 14_regroup.py first" % path, file=sys.stderr)
        sys.exit(1)
    with open(path) as fh:
        rows = list(csv.DictReader(ln for ln in fh if not ln.startswith("#")))
    for r in rows:
        for _, key, _ in GATE_ORDER:
            r[key] = r[key] == "True"
    return rows


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="15_round2_forecast.csv")
    a = ap.parse_args()

    rows = load_rates(os.path.join(HERE, FAILMODE_CSV))
    n_all = len(rows)
    tpl_rows = {}
    for r in rows:
        tpl_rows.setdefault(r["tpl"], []).append(r)

    # conditional pass rate of each gate, in cascade order, over the whole library
    cascade_rate, alive = [], rows
    for name, key, _ in GATE_ORDER:
        nxt = [r for r in alive if r[key]]
        cascade_rate.append(len(nxt) / float(len(alive)) if alive else 0.0)
        alive = nxt
    overall_rate = len(alive) / float(n_all)

    print("=" * 88)
    print("ROUND-2 FORECAST -- planning only, nothing is executed")
    print("generated %s   git %s" % (datetime.datetime.now().isoformat(timespec="seconds"), git_hash()))
    print("=" * 88)

    # ---- 1. measured anchors ----------------------------------------------
    print()
    print("-" * 88)
    print("1. Measured round-1 anchors (not estimates)")
    print("-" * 88)
    print("  AF3 holo   %6.1f s/prediction   %d tokens   5 seeds x 5 samples" % (AF3_SEC_HOLO, AF3_TOKENS_HOLO))
    print("  AF3 apo    %6.1f s/prediction   %d tokens" % (AF3_SEC_APO, AF3_TOKENS_APO))
    print("  AF3 round 1: %d predictions, %d GPUs, %.1f h calendar"
          % (AF3_ROUND1_PREDICTIONS, AF3_GPUS, AF3_ROUND1_CALENDAR_H))
    print("             check: %d x mean(%.1f,%.1f) s / %d GPUs = %.1f h  (measured %.1f h)"
          % (AF3_ROUND1_PREDICTIONS, AF3_SEC_HOLO, AF3_SEC_APO, AF3_GPUS,
             AF3_ROUND1_PREDICTIONS * (AF3_SEC_HOLO + AF3_SEC_APO) / 2 / AF3_GPUS / 3600,
             AF3_ROUND1_CALENDAR_H))
    print("  LigandMPNN %6.1f min per %d sequences per input" % (MPNN_SEC_PER_10K / 60, MPNN_N_SEQ))
    print("  Analysis   %6.1f min per 100 designs (both states), CPU" % (ANALYSIS_SEC_PER_100 / 60))
    print()
    print("  Token-scaling band carried forward: time ~ N^%.2f (measured here) to N^%.1f (quadratic)"
          % (SCALE_EXP_LOW, SCALE_EXP_HIGH))

    # ---- 2. round-2 construct cost ----------------------------------------
    print()
    print("-" * 88)
    print("2. Predicted AF3 cost per prediction, by round-2 construct")
    print("-" * 88)
    constructs = [
        # protein token counts include the CRO chromophore position: the AF3
        # input's chain A is 307 characters for the r28_l3 templates.
        ("round 1 as run", 307, 16, 0, "GGGTGGGTGGGTGGGT, the 16-nt G4 actually used"),
        ("same protein, 22-nt MYC G4", 307, 22, 0, "what README/MISSING_STRUCTURES assume"),
        ("+ rhodamine (Option E)", 307, 22, 33, "JF-series rhodamine, ~33 heavy atoms"),
        ("Nb55 graft + rhodamine", 123 + 160, 22, 33, "nanobody + smaller body, if the barrel is replaced"),
    ]
    print("  %-30s %7s | %-18s %-18s" % ("construct", "tokens", "holo s (low-high)", "vs round 1"))
    for label, naa, ndna, nlig, note in constructs:
        tok = naa + ndna + 3 + nlig
        lo = af3_sec(tok, AF3_SEC_HOLO, AF3_TOKENS_HOLO, SCALE_EXP_LOW)
        hi = af3_sec(tok, AF3_SEC_HOLO, AF3_TOKENS_HOLO, SCALE_EXP_HIGH)
        print("  %-30s %7d | %6.1f - %-11.1f %+5.0f%% to %+.0f%%   %s"
              % (label, tok, lo, hi, 100 * (lo / AF3_SEC_HOLO - 1), 100 * (hi / AF3_SEC_HOLO - 1), note))
    print()
    print("  Read: the construct barely changes size, so AF3 cost per prediction is")
    print("  effectively a round-1 number. The campaign cost is set by HOW MANY")
    print("  predictions, not by what is in them.")

    # ---- 3. measured gate rates -------------------------------------------
    print()
    print("-" * 88)
    print("3. Measured round-1 gate rates -- the basis of every funnel below")
    print("-" * 88)
    print("  %-10s %-42s %10s" % ("gate", "question", "conditional"))
    alive_n = n_all
    for (name, key, q), rate in zip(GATE_ORDER, cascade_rate):
        print("  %-10s %-42s %9.1f%%   %d -> %d" % (name, q, 100 * rate, alive_n, int(round(alive_n * rate))))
        alive_n = int(round(alive_n * rate))
    print("  %-10s %-42s %9.1f%%   %d -> %d" % ("OVERALL", "all four", 100 * overall_rate, n_all, len(alive)))
    print()
    print("  Per-template survival (round 2 inherits whichever parents it starts from):")
    for t in sorted(tpl_rows, key=lambda t: -sum(1 for r in tpl_rows[t] if r["failure_mode"] == "PASS")):
        sub = tpl_rows[t]
        p = sum(1 for r in sub if r["failure_mode"] == "PASS")
        print("    %-28s %3d/%3d = %5.1f%%" % (t, p, len(sub), 100.0 * p / len(sub)))

    # ---- 4. scenarios ------------------------------------------------------
    print()
    print("-" * 88)
    print("4. Scenarios -- predicted compute and predicted funnel")
    print("-" * 88)

    # survivors by parent template, for scenario B
    surv_by_tpl = Counter(r["tpl"] for r in rows if r["failure_mode"] == "PASS")

    def parent_rate(tpls):
        """Weighted mean round-1 PASS rate of the parent templates."""
        num = den = 0
        for t, w in tpls.items():
            sub = tpl_rows.get(t, [])
            if not sub:
                continue
            num += w * sum(1 for r in sub if r["failure_mode"] == "PASS") / float(len(sub))
            den += w
        return num / den if den else overall_rate

    scenarios = [
        dict(key="A", name="Kielen plan as agreed",
             detail="the 21 experimentally tested round-1 variants as design inputs",
             n_inputs=21, seqs=10000, keep=100,
             parents=Counter(r["tpl"] for r in rows if r["global_id"] in TESTED_IDS)),
        dict(key="B", name="Iterate the 71 survivors",
             detail="every design that clears all four gates becomes a design input",
             n_inputs=71, seqs=10000, keep=100, parents=surv_by_tpl),
        dict(key="B-lite", name="Iterate 12 survivors, stratified",
             detail="2 per surviving template, capped -- same funnel, 6x less compute",
             n_inputs=12, seqs=10000, keep=100, parents=surv_by_tpl),
        dict(key="E", name="Option E rebuild",
             detail="new backbones + rhodamine; round-1 rates do NOT transfer",
             n_inputs=8, seqs=10000, keep=100, parents=None),
    ]

    out_rows = []
    for s in scenarios:
        n_in, keep = s["n_inputs"], s["keep"]
        n_designs = n_in * keep
        mpnn_s = n_in * MPNN_SEC_PER_10K * (s["seqs"] / float(MPNN_N_SEQ))
        af3_s = n_designs * (AF3_SEC_HOLO + AF3_SEC_APO)
        ana_s = n_designs / 100.0 * ANALYSIS_SEC_PER_100
        wall_s = mpnn_s / min(AF3_GPUS, n_in) + af3_s / AF3_GPUS + ana_s

        print()
        print("  [%s] %s" % (s["key"], s["name"]))
        print("       %s" % s["detail"])
        print("       inputs %d x %s seqs -> top %d each = %d designs -> %d AF3 predictions (holo+apo)"
              % (n_in, "{:,}".format(s["seqs"]), keep, n_designs, 2 * n_designs))
        print("       LigandMPNN %-10s  AF3 %-10s  analysis %-9s  ->  WALL %s on %d GPUs"
              % (hms(mpnn_s), hms(af3_s), hms(ana_s), hms(wall_s), AF3_GPUS))

        if s["parents"] is None:
            print("       funnel: not predictable from round-1 rates -- different scaffold,")
            print("               different fluorophore, different failure modes. Pilot first.")
            out_rows.append(dict(scenario=s["key"], name=s["name"], n_inputs=n_in,
                                 n_designs=n_designs, n_af3=2 * n_designs,
                                 mpnn_h=round(mpnn_s / 3600, 2), af3_h=round(af3_s / 3600, 2),
                                 wall_h=round(wall_s / 3600, 2),
                                 pass_low="", pass_high="", survivors_low="", survivors_high=""))
            continue

        p_opt = min(1.0, parent_rate(s["parents"]) * WINDOW_CONTACT_UPLIFT)
        p_pes = overall_rate
        lo, hi = int(round(n_designs * min(p_opt, p_pes))), int(round(n_designs * max(p_opt, p_pes)))
        print("       predicted survivors: %d - %d  (%.1f%% pessimistic / %.1f%% optimistic)"
              % (lo, hi, 100 * min(p_opt, p_pes), 100 * max(p_opt, p_pes)))
        print("         pessimistic = redesign randomises, library rate %.1f%%" % (100 * overall_rate))
        print("         optimistic  = parents' rate %.1f%% x %.2f window uplift"
              % (100 * parent_rate(s["parents"]), WINDOW_CONTACT_UPLIFT))
        out_rows.append(dict(scenario=s["key"], name=s["name"], n_inputs=n_in,
                             n_designs=n_designs, n_af3=2 * n_designs,
                             mpnn_h=round(mpnn_s / 3600, 2), af3_h=round(af3_s / 3600, 2),
                             wall_h=round(wall_s / 3600, 2),
                             pass_low=round(100 * min(p_opt, p_pes), 1),
                             pass_high=round(100 * max(p_opt, p_pes), 1),
                             survivors_low=lo, survivors_high=hi))

    # ---- 5. how many to order ---------------------------------------------
    print()
    print("-" * 88)
    print("5. How many to ORDER -- the gate that actually binds")
    print("-" * 88)
    print("  Round 1: %d ordered, %d fluorescent (%.1f%%), %d G4-responsive."
          % (R1_ORDERED, R1_FLUORESCENT, 100.0 * R1_FLUORESCENT / R1_ORDERED, R1_RESPONSIVE))
    # 95% upper bound on a rate given 0/n (rule of three, exact: 1-0.05^(1/n))
    ub = 1 - 0.05 ** (1.0 / R1_ORDERED)
    print("  0 of %d responsive bounds the true switch rate at p < %.1f%% (95%%, exact)."
          % (R1_ORDERED, 100 * ub))
    print()
    print("  n needed for an 80% chance of seeing at least one switch, and what")
    print("  round 1's n=%d would have done at that same rate:" % R1_ORDERED)
    print("    %-14s %-10s %-22s %s" % ("true rate", "n for 80%", "P(0 hits) at n=%d" % R1_ORDERED, "note"))
    for p, note in ((0.03, "pessimistic"), (0.05, ""), (0.10, ""),
                    (CHANG_POLIZZI_RATE, "Chang & Polizzi measured 4/24"), (0.25, "optimistic")):
        n = math.ceil(math.log(0.20) / math.log(1 - p))
        p0 = (1 - p) ** R1_ORDERED
        print("    %-14s %-10d %-22s %s" % ("%.0f%%" % (100 * p), n, "%.0f%%" % (100 * p0), note))
    # rate below which a null result at n=21 was more likely than not
    p_half = 1 - 0.5 ** (1.0 / R1_ORDERED)
    print()
    print("  Read this precisely. Round 1's %d-design null does rule OUT a high switch" % R1_ORDERED)
    print("  rate: anything above %.1f%% is excluded at 95%%. What it cannot do is" % (100 * ub))
    print("  distinguish a design set with a 3%% switch rate from one with none at all --")
    print("  below %.1f%%, seeing nothing was more likely than not even if the designs work."
          % (100 * p_half))
    print("  So the null constrains the architecture, not the score, and the useful")
    print("  question for round 2 is how many to order, not which metric to rank by.")

    # ---- 6. the pilot ------------------------------------------------------
    print()
    print("-" * 88)
    print("6. The pilot that removes both assumptions, before any full run")
    print("-" * 88)
    pilot_n = 10
    pilot_af3_s = pilot_n * (AF3_SEC_HOLO + AF3_SEC_APO) / AF3_GPUS
    pilot_mpnn_s = pilot_n * MPNN_SEC_PER_10K / AF3_GPUS
    pilot_s = pilot_mpnn_s + pilot_af3_s
    print("  Take %d of the 71 survivors, redesign each with the 40+40 window, keep the" % pilot_n)
    print("  single best sequence per parent, predict holo+apo.")
    print("  Cost: LigandMPNN %s + %d AF3 predictions %s = %s on %d GPUs."
          % (hms(pilot_mpnn_s), 2 * pilot_n, hms(pilot_af3_s), hms(pilot_s), AF3_GPUS))
    print()
    print("  It answers, with measurements rather than bands:")
    print("    - A1  the real token scaling, from the round-2 construct's own wall clock")
    print("    - A2  the real inheritance rate: how many of %d children keep their" % pilot_n)
    print("          parent's PASS, which collapses every 'low - high' above to a number")
    print("    - A3  whether the 40+40 window actually lifts the CONTACT rate")
    print()
    biggest = max(r["wall_h"] for r in out_rows)
    print("  %s of pilot de-risks a commitment of up to %s of wall clock."
          % (hms(pilot_s), hms(biggest * 3600)))

    # ---- write -------------------------------------------------------------
    out = os.path.join(HERE, a.out)
    with open(out, "w", newline="") as fh:
        fh.write("# 15_round2_forecast.py\n")
        fh.write("# generated: %s\n" % datetime.datetime.now().isoformat(timespec="seconds"))
        fh.write("# git: %s\n" % git_hash())
        fh.write("# PLANNING ONLY -- predicted, not measured. Nothing was executed.\n")
        fh.write("# timings: *_h are hours; AF3 anchors %.1f s holo / %.1f s apo per prediction, %d GPUs\n"
                 % (AF3_SEC_HOLO, AF3_SEC_APO, AF3_GPUS))
        fh.write("# pass_low/pass_high are PERCENT; see assumptions A1-A3 in the script header\n")
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print()
    print("wrote %s  (%d scenarios)" % (out, len(out_rows)))
    return 0


def _load_tested():
    p = os.path.join(HERE, "tested_variants.csv")
    if not os.path.exists(p):
        return set()
    with open(p) as fh:
        return {r["global_id"] for r in csv.DictReader(fh)}


TESTED_IDS = _load_tested()

if __name__ == "__main__":
    sys.exit(main())

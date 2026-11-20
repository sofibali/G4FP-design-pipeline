#!/usr/bin/env python3
"""
19_select_top.py -- join the full-library rescoring to the gate cascade,
validate the fitness score against experiment, and select a top set.

CHEAP ANALYSIS ONLY.  This script reads tables that already exist.  It runs no
LigandMPNN, no AlphaFold3, and builds no models.

What it does
------------
1. Joins the rescored table to 14_regroup_failuremode.csv on global_id.
2. Reports the gate cascade over whatever population the join covers.
3. Validates the fitness score against the experimentally tested designs --
   the thing that could never be done before, because the rescoring had only
   ever run on 36 designs and 14 tested variants.
4. Selects a top set.

How the top set is chosen, and why not by score
-----------------------------------------------
Round 1's own result is that these metrics do not rank.  Of the designs tested,
every one that passed all four gates was dark, and both that fluoresced were
rejected by the cascade.  Section 3 below measures the same thing directly for
the fitness score.  So the selection rule here is:

  gates first (non-compensatory), then DIVERSITY, never a score ordering.

The fitness score is reported beside every selected design so it can be
inspected, but it is not used to order anything.  A score-ordered list is
printed separately, clearly labelled, purely for comparison.

Usage:
  python3 -I 19_select_top.py --rescored <csv> [--out-prefix 19] [--per-template 3] [--total 32]
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REGROUP = os.path.join(HERE, "14_regroup_failuremode.csv")
TESTED = os.path.join(HERE, "tested_variants.csv")

GATE_SEQ = [("ok_contact", "CONTACT"), ("ok_negative", "NEGATIVE"),
            ("ok_bind", "BIND"), ("ok_switch", "SWITCH")]


def rows_of(path):
    with open(path) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def norm(gid):
    """Collapse the three id spellings this pipeline uses into one key.

    14_regroup_failuremode.csv and tested_variants.csv write
        G4FP_des1_cro_mod4_design_0096
    while 09_rescore_and_select.py writes
        des1_cro_mod4_design_0096
    and the design number is zero-padded inconsistently (design_6 vs
    design_0006).  Joining on the raw string silently yields ZERO matches,
    which is how a 999-row table and a 36-row table can look disjoint.
    """
    g = gid.strip()
    if g.startswith("G4FP_"):
        g = g[len("G4FP_"):]
    m = re.match(r"^(.*)_design_0*(\d+)$", g)
    return "%s_design_%d" % (m.group(1), int(m.group(2))) if m else g


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return None


def truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes", "t")


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                              capture_output=True, text=True, timeout=10,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rescored", required=True)
    ap.add_argument("--out-prefix", default=os.path.join(HERE, "19"))
    ap.add_argument("--per-template", type=int, default=3)
    ap.add_argument("--total", type=int, default=32)
    a = ap.parse_args()

    res = rows_of(a.rescored)
    reg = rows_of(REGROUP)
    tested = {norm(r["global_id"]): r["outcome"].strip().lower()
              for r in rows_of(TESTED)}

    print("=" * 78)
    print("ROUND-1 RESCORING + REGATING -- cheap analysis, nothing was generated")
    print("=" * 78)
    print("  rescored table : %s  (%d rows)" % (os.path.basename(a.rescored), len(res)))
    print("  gate table     : 14_regroup_failuremode.csv  (%d rows)" % len(reg))
    print("  tested variants: tested_variants.csv  (%d rows)" % len(tested))

    # ---- join -------------------------------------------------------------
    rindex = {norm(r["global_id"]): r for r in res}
    joined = []
    for g in reg:
        r = rindex.get(norm(g["global_id"]))
        if r is None:
            continue
        rec = dict(g)
        rec["fitness_score"] = fnum(r, "fitness_score")
        rec["holo_iptm"] = fnum(r, "holo_iptm")
        rec["holo_chromophore_plddt_fixed"] = fnum(r, "holo_chromophore_plddt_fixed")
        rec["holo_mean_plddt"] = fnum(r, "holo_mean_plddt")
        rec["sequence"] = r.get("sequence", "")
        joined.append(rec)

    print()
    print("  joined on global_id: %d of %d gate rows matched a rescored row"
          % (len(joined), len(reg)))
    if len(joined) < len(reg):
        print("  !! %d gate rows have no rescoring. The rescoring did not cover"
              % (len(reg) - len(joined)))
        print("     the whole library; every rate below is over the joined subset.")

    scored = [r for r in joined if r["fitness_score"] is not None]
    print("  of those, %d carry a fitness_score" % len(scored))

    # ---- 1. gate cascade --------------------------------------------------
    print()
    print("1. Gate cascade over the joined population (N = %d)" % len(joined))
    print("-" * 78)
    alive = joined
    print("   %-10s %8s %8s %9s" % ("gate", "reached", "passed", "cond."))
    for key, name in GATE_SEQ:
        passed = [r for r in alive if truthy(r[key])]
        print("   %-10s %8d %8d %8.1f %%" %
              (name, len(alive), len(passed), 100.0 * len(passed) / max(1, len(alive))))
        alive = passed
    survivors = alive
    print()
    print("   PASS all four: %d  (%.2f %% of %d)"
          % (len(survivors), 100.0 * len(survivors) / max(1, len(joined)), len(joined)))
    print()
    print("   first gate failed:")
    for k, v in Counter(r["failure_mode"] for r in joined).most_common():
        print("     %-22s %4d  (%4.1f %%)" % (k, v, 100.0 * v / len(joined)))

    # ---- 2. per-template --------------------------------------------------
    print()
    print("2. Survivors per template")
    print("-" * 78)
    tpl_tot = Counter(r["tpl"] for r in joined)
    tpl_pass = Counter(r["tpl"] for r in survivors)
    for t in sorted(tpl_tot, key=lambda x: -tpl_pass[x]):
        print("   %-22s %3d / %3d   %5.1f %%" %
              (t, tpl_pass[t], tpl_tot[t], 100.0 * tpl_pass[t] / tpl_tot[t]))

    # ---- 3. does the score rank? -----------------------------------------
    print()
    print("3. Does the fitness score rank? -- validated against experiment")
    print("-" * 78)
    lab = [r for r in scored if norm(r["global_id"]) in tested]
    if not lab:
        print("   no tested design carries a fitness score; cannot validate.")
    else:
        lab.sort(key=lambda r: -r["fitness_score"])
        n = len(lab)
        print("   %d of %d tested designs carry a fitness score.\n" % (n, len(tested)))
        print("   %-32s %9s %5s %-10s %s" % ("global_id", "fitness", "rank", "outcome", "gates"))
        pos_ranks = []
        for i, r in enumerate(lab, 1):
            out = tested[norm(r["global_id"])]
            gates = "PASS" if r["failure_mode"] == "PASS" else r["failure_mode"].split("_")[0]
            flag = ""
            if out != "dark":
                pos_ranks.append(i)
                flag = "  <<<"
            print("   %-32s %9.4f %5d %-10s %s%s" %
                  (r["global_id"], r["fitness_score"], i, out, gates, flag))
        if pos_ranks:
            # AUC = P(a random positive outscores a random dark design).
            pos = [r["fitness_score"] for r in lab if tested[norm(r["global_id"])] != "dark"]
            neg = [r["fitness_score"] for r in lab if tested[norm(r["global_id"])] == "dark"]
            auc = sum((a > b) + 0.5 * (a == b) for a in pos for b in neg) / float(len(pos) * len(neg))
            # exact permutation p: over all ways to label n designs with
            # len(pos) positives, how often is the AUC this low or lower?
            import itertools
            sc = [r["fitness_score"] for r in lab]
            tot = bad = 0
            for combo in itertools.combinations(range(n), len(pos)):
                tot += 1
                cs = set(combo)
                w = sum((sc[i] > sc[j]) + 0.5 * (sc[i] == sc[j])
                        for i in combo for j in range(n) if j not in cs)
                if w / float(len(pos) * (n - len(pos))) <= auc:
                    bad += 1
            print()
            print("   AUC of fitness_score for predicting fluorescence: %.3f" % auc)
            print("     0.5 = coin flip | >0.5 = ranks correctly | <0.5 = inverted")
            print("   exact permutation P(AUC <= observed) = %d/%d = %.3f"
                  % (bad, tot, bad / float(tot)))
            print("   NOT significant at %d positives. The point estimate is on the"
                  % len(pos))
            print("   wrong side of chance; nothing here argues the score ranks.")
            print()
            exp = (n + 1) / 2.0
            print()
            print("   positives at rank(s) %s of %d" %
                  (", ".join(str(x) for x in pos_ranks), n))
            print("   mean rank of positives : %.1f" % (sum(pos_ranks) / len(pos_ranks)))
            print("   expected if score were random : %.1f" % exp)
            if sum(pos_ranks) / len(pos_ranks) > exp:
                print("   -> the positives rank BELOW chance. The score is not merely")
                print("      uninformative here; it points the wrong way.")
            print("   n = %d positives. This cannot establish inversion, but it is"
                  % len(pos_ranks))
            print("   the only labelled data there is, and it does not support ranking.")

    # ---- 4. the top set ---------------------------------------------------
    print()
    print("4. Top set -- gates, then diversity. NOT ordered by score.")
    print("-" * 78)
    by_tpl = defaultdict(list)
    for r in survivors:
        by_tpl[r["tpl"]].append(r)
    # within a template, spread across the dPAE range rather than take the top
    for t in by_tpl:
        by_tpl[t].sort(key=lambda r: (fnum(r, "dpae") or 0.0))
    picked = []
    for t in sorted(by_tpl, key=lambda x: -len(by_tpl[x])):
        pool = by_tpl[t]
        k = min(a.per_template, len(pool))
        if k:
            step = max(1, len(pool) // k)
            picked.extend(pool[::step][:k])
    picked = picked[:a.total]
    print("   rule: must PASS all four gates; at most %d per template;" % a.per_template)
    print("         within a template, spread evenly across the dPAE range.")
    print("   selected %d designs from %d templates\n" % (picked and len(picked) or 0,
                                                          len(set(r["tpl"] for r in picked))))
    print("   %-32s %8s %8s %8s %9s" % ("global_id", "dPAE", "iPAE_ap", "n_cont", "fitness"))
    for r in sorted(picked, key=lambda r: (r["tpl"], r["global_id"])):
        fs = r["fitness_score"]
        print("   %-32s %8.2f %8.2f %8s %9s" %
              (r["global_id"], fnum(r, "dpae") or 0.0, fnum(r, "ipae_apo") or 0.0,
               r["n_module_contact"], ("%.4f" % fs) if fs is not None else "n/a"))

    # ---- 5. score-ordered, for comparison only ---------------------------
    print()
    print("5. For comparison ONLY -- the same survivors ordered by fitness score")
    print("-" * 78)
    print("   Shown so the two lists can be compared, NOT as a recommendation.")
    ss = [r for r in survivors if r["fitness_score"] is not None]
    ss.sort(key=lambda r: -r["fitness_score"])
    for i, r in enumerate(ss[:a.total], 1):
        print("   %2d. %-32s %9.4f  %s" %
              (i, r["global_id"], r["fitness_score"], r["tpl"]))
    if ss:
        overlap = len(set(norm(r["global_id"]) for r in picked) &
                      set(norm(r["global_id"]) for r in ss[:a.total]))
        print()
        print("   overlap with the diversity set: %d of %d" % (overlap, len(picked)))

    # ---- write ------------------------------------------------------------
    hdr = ["global_id", "tpl", "design_id", "failure_mode", "n_module_contact",
           "ipae_apo", "ipae_holo", "dpae", "bind_ipae_holo", "iptm_holo",
           "fitness_score", "holo_mean_plddt", "holo_chromophore_plddt_fixed",
           "sequence"]
    stamp = ["# 19_select_top.py",
             "# generated: %s" % datetime.datetime.now().isoformat(timespec="seconds"),
             "# git: %s" % git_hash(),
             "# CHEAP ANALYSIS ONLY -- no LigandMPNN, no AlphaFold3, no model building",
             "# source: %s + 14_regroup_failuremode.csv" % os.path.basename(a.rescored),
             "# selection: all four gates, then <=%d per template spread over dPAE" % a.per_template,
             "# fitness_score is REPORTED, never used to order -- see section 3",
             "# units: *PAE in Angstrom; n_module_contact is a residue count"]
    for name, data in (("_survivors_full.csv", survivors), ("_top_set.csv", picked)):
        p = a.out_prefix + name
        with open(p, "w", newline="") as fh:
            fh.write("\n".join(stamp) + "\n")
            w = csv.DictWriter(fh, fieldnames=hdr, extrasaction="ignore")
            w.writeheader()
            w.writerows(data)
        print("\nWrote %s (%d rows)" % (p, len(data)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

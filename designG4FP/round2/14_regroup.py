#!/usr/bin/env python3
"""
14_regroup.py -- regroup all 999 round-1 designs by FIRST GATE FAILED.

The cascade in 13_RESCORING_PLAN.md reports how many survive each stage.  That
is not the same question as "why did each design die".  A design that fails
CONTACT is thrown away before NEGATIVE, BIND or SWITCH ever see it, so the
stage-by-stage counts tell you nothing about how the 701 non-survivors are
distributed across causes -- and the cause is what round 2 has to fix.

This script assigns every design exactly one failure mode (the first gate it
fails, in cascade order) plus the independent "would it have passed gate X on
its own" flags, and regroups by template.

Gate order and thresholds are read from, and must match, 13_RESCORING_PLAN.md.

  CONTACT   n_module_contact >= 3        does AF3 think the module touches the G4
  NEGATIVE  ipae_apo        >= 8.0 A     is there anything left to close
  BIND      bind_ipae_holo  <= P40       is the protein-G4 interface credible
  SWITCH    dpae            >  3.0 A     Chang & Polizzi's published threshold

No re-prediction, no new structures: this reads 13_delta_ipae_module.csv, which
was itself computed from PAE matrices AlphaFold3 wrote in March 2026.

Usage:  python3 -I 14_regroup.py [--gate-csv 13_delta_ipae_module.csv]
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CONTACT_MIN = 3
APO_IPAE_MIN = 8.0
BIND_PCTILE = 40.0
SWITCH_MIN = 3.0

FLOAT_COLS = ("n_module_contact", "ipae_apo", "ipae_holo", "dpae",
              "bind_ipae_holo", "bind_pae_min", "iptm_holo", "ptm_holo")
# dpae_contact is deliberately NOT required: it is empty by construction for the
# 152 designs with n_module_contact == 0, because there are no contact residues
# to average over.  Dropping those rows would silently delete the worst-scoring
# sixth of the library and inflate every pass rate below.
OPTIONAL_FLOAT_COLS = ("dpae_contact",)


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                              capture_output=True, text=True, timeout=10).stdout.strip() or "n/a"
    except Exception:
        return "n/a"


def read_gate(path):
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows, dropped = [], []
    for r in csv.DictReader(lines):
        tpl = r["template"][len("output_"):] if r["template"].startswith("output_") else r["template"]
        r["tpl"] = tpl
        r["global_id"] = "%s_design_%04d" % (tpl, int(r["design_id"].rsplit("_", 1)[1]))
        ok = True
        for k in FLOAT_COLS:
            try:
                r[k] = float(r[k])
            except (TypeError, ValueError):
                r[k] = float("nan")
                ok = False
        for k in OPTIONAL_FLOAT_COLS:
            try:
                r[k] = float(r[k])
            except (TypeError, ValueError):
                r[k] = float("nan")
        if ok:
            rows.append(r)
        else:
            dropped.append(r["global_id"])
    if dropped:
        print("NOTE: dropped %d rows with unparseable required values: %s"
              % (len(dropped), ", ".join(dropped[:5])), file=sys.stderr)
    return rows


def read_tested(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["global_id"]: r["outcome"] for r in csv.DictReader(fh)}


def percentile(sorted_vals, p):
    if not sorted_vals:
        return float("nan")
    i = min(len(sorted_vals) - 1, max(0, int(round(p / 100.0 * (len(sorted_vals) - 1)))))
    return sorted_vals[i]


def classify(rows, bind_cut):
    """First gate failed, in cascade order.  PASS if it clears all four."""
    for r in rows:
        if r["n_module_contact"] < CONTACT_MIN:
            r["failure_mode"] = "CONTACT_no_interface"
        elif r["ipae_apo"] < APO_IPAE_MIN:
            r["failure_mode"] = "NEGATIVE_pre_closed"
        elif r["bind_ipae_holo"] > bind_cut:
            r["failure_mode"] = "BIND_weak_interface"
        elif r["dpae"] <= SWITCH_MIN:
            r["failure_mode"] = "SWITCH_no_contrast"
        else:
            r["failure_mode"] = "PASS"
        # independent flags: would it pass this gate on its own?
        r["ok_contact"] = r["n_module_contact"] >= CONTACT_MIN
        r["ok_negative"] = r["ipae_apo"] >= APO_IPAE_MIN
        r["ok_bind"] = r["bind_ipae_holo"] <= bind_cut
        r["ok_switch"] = r["dpae"] > SWITCH_MIN
    return rows


MODES = ["CONTACT_no_interface", "NEGATIVE_pre_closed", "BIND_weak_interface",
         "SWITCH_no_contrast", "PASS"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate-csv", default="13_delta_ipae_module.csv")
    ap.add_argument("--out-prefix", default="14_regroup")
    a = ap.parse_args()

    gate_path = a.gate_csv if os.path.isabs(a.gate_csv) else os.path.join(HERE, a.gate_csv)
    rows = read_gate(gate_path)
    tested = read_tested(os.path.join(HERE, "tested_variants.csv"))

    # 13_RESCORING_PLAN.md section 4: thresholds are percentiles of the measured
    # ROUND-1 DISTRIBUTION, i.e. of the whole library, not of the subset that
    # happens to reach this gate.  That keeps the threshold a fixed property of
    # the population, so it does not move when an earlier gate is retuned.
    bind_cut = percentile(sorted(r["bind_ipae_holo"] for r in rows), BIND_PCTILE)
    reach_bind = [r for r in rows if r["n_module_contact"] >= CONTACT_MIN
                  and r["ipae_apo"] >= APO_IPAE_MIN]
    bind_cut_conditional = percentile(sorted(r["bind_ipae_holo"] for r in reach_bind), BIND_PCTILE)

    classify(rows, bind_cut)

    print("=" * 86)
    print("REGROUPING -- %d designs, first gate failed" % len(rows))
    print("=" * 86)
    print("  gate csv     : %s" % os.path.basename(gate_path))
    print("  BIND cut     : interface PAE <= %.2f A  (P%.0f of all %d, per plan section 4)"
          % (bind_cut, BIND_PCTILE, len(rows)))
    print("                 for comparison, P%.0f of only the %d that reach this gate "
          "would be %.2f A" % (BIND_PCTILE, len(reach_bind), bind_cut_conditional))
    print()

    # ---- 1. failure-mode census -------------------------------------------
    print("-" * 86)
    print("1. Where the library dies")
    print("-" * 86)
    print("  %-24s %6s %8s   %s" % ("first gate failed", "n", "share", "what it means"))
    meaning = {
        "CONTACT_no_interface": "AF3 predicts <3 module-G4 contacts",
        "NEGATIVE_pre_closed": "module already shut with no G4 present",
        "BIND_weak_interface": "interface PAE in the worst 60%",
        "SWITCH_no_contrast": "no apo/holo confidence contrast",
        "PASS": "clears all four",
    }
    for m in MODES:
        n = sum(1 for r in rows if r["failure_mode"] == m)
        print("  %-24s %6d %7.1f%%   %s" % (m, n, 100.0 * n / len(rows), meaning[m]))

    # ---- 2. independent pass rates ----------------------------------------
    print()
    print("-" * 86)
    print("2. Each gate ALONE, over all %d (not conditioned on the ones before it)" % len(rows))
    print("-" * 86)
    for key, label in (("ok_contact", "CONTACT  >=3 contacts"),
                       ("ok_negative", "NEGATIVE apo iPAE >= 8 A"),
                       ("ok_bind", "BIND     iPAE <= P40"),
                       ("ok_switch", "SWITCH   dPAE > 3 A")):
        n = sum(1 for r in rows if r[key])
        print("  %-28s %4d / %d  = %5.1f%%" % (label, n, len(rows), 100.0 * n / len(rows)))
    both = sum(1 for r in rows if r["ok_switch"] and r["ok_bind"])
    print("  %-28s %4d / %d  = %5.1f%%  (expected %.0f if independent)"
          % ("SWITCH and BIND together", both, len(rows), 100.0 * both / len(rows),
             sum(1 for r in rows if r["ok_switch"]) * sum(1 for r in rows if r["ok_bind"]) / len(rows)))

    # ---- 3. regroup by template -------------------------------------------
    print()
    print("-" * 86)
    print("3. By template -- the cascade is not uniform across the ten inputs")
    print("-" * 86)
    tpls = sorted({r["tpl"] for r in rows})
    hdr = "  %-26s %4s | %5s %5s %5s %5s | %5s %6s" % (
        "template", "n", "CONT", "PRECL", "WEAK", "NOSW", "PASS", "rate")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for t in tpls:
        sub = [r for r in rows if r["tpl"] == t]
        c = {m: sum(1 for r in sub if r["failure_mode"] == m) for m in MODES}
        print("  %-26s %4d | %5d %5d %5d %5d | %5d %5.1f%%"
              % (t, len(sub), c["CONTACT_no_interface"], c["NEGATIVE_pre_closed"],
                 c["BIND_weak_interface"], c["SWITCH_no_contrast"], c["PASS"],
                 100.0 * c["PASS"] / len(sub)))
    c = {m: sum(1 for r in rows if r["failure_mode"] == m) for m in MODES}
    print("  " + "-" * (len(hdr) - 2))
    print("  %-26s %4d | %5d %5d %5d %5d | %5d %5.1f%%"
          % ("ALL", len(rows), c["CONTACT_no_interface"], c["NEGATIVE_pre_closed"],
             c["BIND_weak_interface"], c["SWITCH_no_contrast"], c["PASS"],
             100.0 * c["PASS"] / len(rows)))

    # ---- 4. per-template medians, so round 2 can pick inputs ---------------
    print()
    print("-" * 86)
    print("4. Per-template medians -- which scaffolds are worth redesigning")
    print("-" * 86)
    print("  %-26s %8s %8s %8s %8s %8s" % ("template", "contacts", "apo iPAE", "dPAE", "bind iPAE", "iPTM"))
    for t in tpls:
        sub = [r for r in rows if r["tpl"] == t]
        med = lambda k: statistics.median([r[k] for r in sub])
        print("  %-26s %8.1f %8.2f %8.2f %8.2f %8.3f"
              % (t, med("n_module_contact"), med("ipae_apo"), med("dpae"),
                 med("bind_ipae_holo"), med("iptm_holo")))

    # ---- 5. the tested 21, by failure mode --------------------------------
    if tested:
        print()
        print("-" * 86)
        print("5. The 21 experimentally tested designs, regrouped")
        print("-" * 86)
        by_id = {r["global_id"]: r for r in rows}
        found = [(g, o, by_id[g]) for g, o in tested.items() if g in by_id]
        for m in MODES:
            sub = [(g, o) for g, o, r in found if r["failure_mode"] == m]
            if not sub:
                continue
            nb = sum(1 for _, o in sub if o.startswith("fluorescent"))
            print("  %-24s %2d tested  (%d fluorescent)" % (m, len(sub), nb))
            for g, o in sorted(sub):
                print("      %-36s %s" % (g, o))

    # ---- 6. survivors ------------------------------------------------------
    surv = [r for r in rows if r["failure_mode"] == "PASS"]
    surv.sort(key=lambda r: -r["dpae"])
    print()
    print("-" * 86)
    print("6. The %d survivors, by template" % len(surv))
    print("-" * 86)
    from collections import Counter
    for t, n in Counter(r["tpl"] for r in surv).most_common():
        print("  %-26s %3d" % (t, n))
    print()
    print("  top 10 by dPAE:")
    print("  %-36s %7s %9s %9s %6s" % ("global_id", "dPAE", "bind iPAE", "apo iPAE", "cont"))
    for r in surv[:10]:
        print("  %-36s %7.2f %9.2f %9.2f %6.0f"
              % (r["global_id"], r["dpae"], r["bind_ipae_holo"], r["ipae_apo"], r["n_module_contact"]))

    # ---- write -------------------------------------------------------------
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    head = ["# 14_regroup.py", "# generated: %s" % stamp, "# git: %s" % git_hash(),
            "# source: %s" % os.path.basename(gate_path),
            "# gates: CONTACT n>=%d | NEGATIVE apo iPAE>=%.1f | BIND iPAE<=%.2f (P%.0f) | SWITCH dPAE>%.1f"
            % (CONTACT_MIN, APO_IPAE_MIN, bind_cut, BIND_PCTILE, SWITCH_MIN),
            "# failure_mode = FIRST gate failed, in that order",
            "# units: all *PAE columns in Angstrom; n_module_contact is a residue count"]
    cols = ["global_id", "tpl", "design_id", "failure_mode", "n_module_contact", "ipae_apo",
            "ipae_holo", "dpae", "dpae_contact", "bind_ipae_holo", "bind_pae_min",
            "iptm_holo", "ptm_holo", "ok_contact", "ok_negative", "ok_bind", "ok_switch"]

    def dump(path, data):
        with open(path, "w", newline="") as fh:
            for h in head:
                fh.write(h + "\n")
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in data:
                w.writerow(r)
        print("wrote %s  (%d rows)" % (path, len(data)))

    print()
    dump(os.path.join(HERE, a.out_prefix + "_failuremode.csv"), rows)
    dump(os.path.join(HERE, a.out_prefix + "_survivors.csv"), surv)
    return 0


if __name__ == "__main__":
    sys.exit(main())

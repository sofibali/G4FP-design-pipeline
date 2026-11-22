#!/usr/bin/env python3
"""Diff the committed (pre-fix) result tables against the corrected re-run.

The point of this script is that the *delta* is the result, not the corrected
table.  A corrected table published on its own silently replaces the old
numbers; anyone who already acted on them has no way to see what moved.  So
before committing the re-run we report, for every changed table:

  * the mean shift of each numeric column, old -> new;
  * the Pearson r between the old and new values of that column (a column
    that moved by a constant offset still has r = 1; a column that was
    recomputed from scratch usually does not);
  * for the master scored table, the rank stability of `fitness_score`
    (Spearman rho, top-10 and top-50 overlap), because a ranking is what the
    table is actually used for;
  * where the two experimentally fluorescent designs moved.

"Old" is read straight out of git (HEAD:<path>), so this script requires the
re-run to be uncommitted -- which is exactly when you need it.  Stdlib only;
the statistics here are simple enough that an extra dependency would be a
bigger risk than the arithmetic.

Usage:
    python3 -I 20_delta_rerun.py --repo /path/to/repo --out 20_RERUN_DELTA.md
"""

import argparse
import csv
import io
import math
import os
import subprocess
import sys
from datetime import datetime, timezone

# Columns that are identifiers or flags, not measurements.  Shifts in these
# are reported as set changes, not as means.
ID_COLS = {"seq_id", "design_name", "template", "state", "sequence",
           "global_id", "is_template"}
FLAG_COLS = {"in_pareto", "in_top10_per_tpl", "in_top50_global"}

MASTER = "designG4FP/07_results/plot_data_csvs/07_all_designs_scored.csv"


def git(repo, *args):
    return subprocess.run(["git", "-C", repo] + list(args),
                          check=True, capture_output=True, text=True).stdout


def read_csv_text(text):
    return list(csv.DictReader(io.StringIO(text)))


def changed_csvs(repo):
    out = git(repo, "status", "--porcelain")
    paths = []
    for line in out.splitlines():
        p = line[3:].strip()
        if p.endswith(".csv"):
            paths.append(p)
    return sorted(paths)


def to_float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def key_of(row):
    """Join key: global_id when the table has one, else seq_id (+state)."""
    if row.get("global_id"):
        return row["global_id"]
    k = row.get("seq_id") or row.get("design_name") or ""
    if row.get("state"):
        k = "%s|%s" % (k, row["state"])
    return k


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return None          # one side is constant; r is undefined, not 1
    return sxy / math.sqrt(sxx * syy)


def ranks(vals):
    """Tie-averaged ranks, 1 = largest (these are scores: bigger is better)."""
    order = sorted(range(len(vals)), key=lambda i: -vals[i])
    out = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def spearman(xs, ys):
    return pearson(ranks(xs), ranks(ys))


def compare_table(path, old_rows, new_rows):
    """Per-column delta for one table.  Returns a dict of findings."""
    old = {key_of(r): r for r in old_rows}
    new = {key_of(r): r for r in new_rows}
    shared = [k for k in new if k in old]

    cols = [c for c in (new_rows[0].keys() if new_rows else [])
            if c not in ID_COLS]
    stats = []
    for c in cols:
        pairs = []
        for k in shared:
            a, b = to_float(old[k].get(c)), to_float(new[k].get(c))
            if a is not None and b is not None:
                pairs.append((a, b))
        if len(pairs) < 3:
            continue
        xs = [a for a, _ in pairs]
        ys = [b for _, b in pairs]
        n_changed = sum(1 for a, b in pairs if abs(a - b) > 1e-9)
        if n_changed == 0:
            continue
        mo = sum(xs) / len(xs)
        mn = sum(ys) / len(ys)
        stats.append({
            "column": c,
            "n": len(pairs),
            "n_changed": n_changed,
            "mean_old": mo,
            "mean_new": mn,
            "shift": mn - mo,
            "pct": (100.0 * (mn - mo) / abs(mo)) if abs(mo) > 1e-12 else None,
            "r": pearson(xs, ys),
            "max_abs_change": max(abs(a - b) for a, b in pairs),
        })
    stats.sort(key=lambda s: -abs(s["shift"]))
    return {
        "path": path,
        "n_old": len(old_rows),
        "n_new": len(new_rows),
        "n_shared": len(shared),
        "only_old": sorted(set(old) - set(new)),
        "only_new": sorted(set(new) - set(old)),
        "stats": stats,
    }


def master_ranking(old_rows, new_rows, positives):
    """Rank stability of fitness_score, and where the positives moved."""
    old = {key_of(r): r for r in old_rows}
    new = {key_of(r): r for r in new_rows}
    shared = sorted(set(old) & set(new))
    xs, ys, keys = [], [], []
    for k in shared:
        a, b = to_float(old[k].get("fitness_score")), to_float(new[k].get("fitness_score"))
        if a is not None and b is not None:
            keys.append(k); xs.append(a); ys.append(b)

    ro = dict(zip(keys, ranks(xs)))
    rn = dict(zip(keys, ranks(ys)))

    def topn(rmap, n):
        return set(sorted(rmap, key=lambda k: rmap[k])[:n])

    res = {
        "n": len(keys),
        "pearson": pearson(xs, ys),
        "spearman": spearman(xs, ys),
        "top10_overlap": len(topn(ro, 10) & topn(rn, 10)),
        "top50_overlap": len(topn(ro, 50) & topn(rn, 50)),
        "top100_overlap": len(topn(ro, 100) & topn(rn, 100)),
        "positives": [],
        "biggest_movers": [],
    }
    for p in positives:
        if p in ro:
            res["positives"].append({
                "global_id": p,
                "rank_old": int(round(ro[p])), "rank_new": int(round(rn[p])),
                "score_old": to_float(old[p].get("fitness_score")),
                "score_new": to_float(new[p].get("fitness_score")),
            })
    movers = sorted(keys, key=lambda k: -abs(rn[k] - ro[k]))[:10]
    for k in movers:
        res["biggest_movers"].append({
            "global_id": k,
            "rank_old": int(round(ro[k])), "rank_new": int(round(rn[k])),
            "delta": int(round(rn[k] - ro[k])),
        })
    return res


def per_template_change(old_rows, new_rows, columns):
    """How many designs changed, per template, for the named columns.

    A correction that lands on some templates and not others is a different
    kind of finding from one that shifts everything a little: it says the bug
    was conditional on something template-specific (residue numbering, chain
    order), and it tells you which earlier conclusions are untouched.
    """
    old = {key_of(r): r for r in old_rows}
    out = {}
    for c in columns:
        tally = {}
        for r in new_rows:
            k = key_of(r)
            if k not in old:
                continue
            t = r.get("template", "?")
            a, b = to_float(old[k].get(c)), to_float(r.get(c))
            if a is None or b is None:
                continue
            d = tally.setdefault(t, [0, 0])
            d[1] += 1
            if abs(a - b) > 1e-9:
                d[0] += 1
        out[c] = tally
    return out


def flag_changes(old_rows, new_rows):
    old = {key_of(r): r for r in old_rows}
    new = {key_of(r): r for r in new_rows}
    shared = set(old) & set(new)
    out = {}
    for c in FLAG_COLS:
        if not new_rows or c not in new_rows[0]:
            continue
        truthy = lambda v: str(v).strip().lower() in ("1", "true", "yes")
        so = {k for k in shared if truthy(old[k].get(c))}
        sn = {k for k in shared if truthy(new[k].get(c))}
        out[c] = {"n_old": len(so), "n_new": len(sn),
                  "kept": len(so & sn),
                  "dropped": sorted(so - sn), "added": sorted(sn - so)}
    return out


def fmt(v, nd=4):
    if v is None:
        return "n/a"
    return ("%%.%df" % nd) % v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--tested", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    repo = os.path.abspath(a.repo)
    head = git(repo, "rev-parse", "HEAD").strip()
    paths = changed_csvs(repo)
    if not paths:
        sys.exit("No modified CSVs in the working tree -- nothing to diff. "
                 "This script must run BEFORE the re-run is committed.")

    positives = []
    with open(a.tested) as fh:
        for r in csv.DictReader(fh):
            if r.get("outcome", "").strip().lower() != "dark":
                positives.append(r["global_id"].strip())

    reports = []
    master_rep = None
    for p in paths:
        try:
            old_rows = read_csv_text(git(repo, "show", "HEAD:%s" % p))
        except subprocess.CalledProcessError:
            continue                      # new file, not a re-run
        with open(os.path.join(repo, p)) as fh:
            new_rows = list(csv.DictReader(fh))
        if not old_rows or not new_rows:
            continue
        rep = compare_table(p, old_rows, new_rows)
        reports.append(rep)
        if p == MASTER:
            master_rep = (rep, master_ranking(old_rows, new_rows, positives),
                          flag_changes(old_rows, new_rows),
                          per_template_change(old_rows, new_rows,
                              ["chromophore_rmsd", "holo_chromophore_plddt",
                               "fitness_score"]))

    L = []
    w = L.append
    w("# Re-run delta: what the corrected analysis changed")
    w("")
    w("Generated %s by `round2/20_delta_rerun.py`."
      % datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    w("Baseline is git `%s` (the committed, pre-fix tables); \"new\" is the"
      % head[:9])
    w("working tree after the corrected re-run. %d changed CSV tables."
      % len(reports))
    w("")
    w("**Read the shifts, not just the new tables.** A corrected table that")
    w("replaces the old one in place hides the size of the correction from")
    w("anyone who already used the old numbers.")
    w("")

    if master_rep:
        rep, rk, flags, ptc = master_rep
        w("## 1. The master scored table")
        w("")
        w("`%s`" % MASTER)
        w("")
        w("%d rows before, %d after, %d joined on `global_id`."
          % (rep["n_old"], rep["n_new"], rep["n_shared"]))
        w("")
        w("### 1.1 Ranking stability of `fitness_score`")
        w("")
        w("| quantity | value |")
        w("|---|---:|")
        w("| designs compared | %d |" % rk["n"])
        w("| Pearson r (old vs new score) | %s |" % fmt(rk["pearson"], 3))
        w("| Spearman rho (old vs new ordering) | %s |" % fmt(rk["spearman"], 3))
        w("| top-10 overlap | %d / 10 |" % rk["top10_overlap"])
        w("| top-50 overlap | %d / 50 |" % rk["top50_overlap"])
        w("| top-100 overlap | %d / 100 |" % rk["top100_overlap"])
        w("")
        if rk["positives"]:
            w("### 1.2 Where the two fluorescent designs moved")
            w("")
            w("| design | rank before | rank after | score before | score after |")
            w("|---|---:|---:|---:|---:|")
            for p in rk["positives"]:
                w("| `%s` | %d | %d | %s | %s |"
                  % (p["global_id"], p["rank_old"], p["rank_new"],
                     fmt(p["score_old"], 4), fmt(p["score_new"], 4)))
            w("")
            w("Ranks are out of %d. These are the only two designs in the"
              % rk["n"])
            w("campaign with a measured positive phenotype, so where they sit")
            w("is the only part of the ranking that has an external check.")
            w("")
        w("### 1.3 Selection-flag churn")
        w("")
        w("| flag | before | after | kept | dropped | added |")
        w("|---|---:|---:|---:|---:|---:|")
        for c, f in sorted(flags.items()):
            w("| `%s` | %d | %d | %d | %d | %d |"
              % (c, f["n_old"], f["n_new"], f["kept"],
                 len(f["dropped"]), len(f["added"])))
        w("")
        w("### 1.5 Which templates the correction touched")
        w("")
        cols = list(ptc.keys())
        w("| template | " + " | ".join("`%s`" % c for c in cols) + " |")
        w("|---" * (len(cols) + 1) + "|")
        tnames = sorted({t for c in cols for t in ptc[c]})
        for t in tnames:
            cells = []
            for c in cols:
                ch, tot = ptc[c].get(t, (0, 0))
                cells.append("%d / %d" % (ch, tot))
            w("| `%s` | %s |" % (t, " | ".join(cells)))
        w("")
        w("A column that changed on some templates and not others means the")
        w("defect was conditional on something template-specific, so the")
        w("untouched templates' earlier numbers still stand.")
        w("")
        w("### 1.4 Columns that moved")
        w("")
        w("| column | n changed | mean before | mean after | shift | r(old,new) | max abs change |")
        w("|---|---:|---:|---:|---:|---:|---:|")
        for s in rep["stats"]:
            w("| `%s` | %d | %s | %s | %s | %s | %s |"
              % (s["column"], s["n_changed"], fmt(s["mean_old"], 3),
                 fmt(s["mean_new"], 3), fmt(s["shift"], 3),
                 fmt(s["r"], 3), fmt(s["max_abs_change"], 3)))
        w("")

    w("## 2. Every changed table")
    w("")
    w("One line per table; the column with the largest absolute mean shift.")
    w("")
    w("| table | rows | cols changed | largest shift | mean before | mean after | r |")
    w("|---|---:|---:|---|---:|---:|---:|")
    for rep in sorted(reports, key=lambda r: r["path"]):
        if rep["stats"]:
            s = rep["stats"][0]
            w("| `%s` | %d | %d | `%s` %s | %s | %s | %s |"
              % (rep["path"], rep["n_new"], len(rep["stats"]), s["column"],
                 fmt(s["shift"], 3), fmt(s["mean_old"], 3),
                 fmt(s["mean_new"], 3), fmt(s["r"], 3)))
        else:
            w("| `%s` | %d | 0 | (no numeric change) | | | |"
              % (rep["path"], rep["n_new"]))
    w("")

    moved = [r for r in reports if r["only_old"] or r["only_new"]]
    if moved:
        w("## 3. Tables whose row membership changed")
        w("")
        for rep in moved:
            w("- `%s`: %d rows dropped, %d added"
              % (rep["path"], len(rep["only_old"]), len(rep["only_new"])))
        w("")

    if master_rep:
        w("## 4. Biggest individual rank moves")
        w("")
        w("| design | rank before | rank after | delta |")
        w("|---|---:|---:|---:|")
        for m in master_rep[1]["biggest_movers"]:
            w("| `%s` | %d | %d | %+d |"
              % (m["global_id"], m["rank_old"], m["rank_new"], m["delta"]))
        w("")

    w("## 5. How to reproduce this")
    w("")
    w("```")
    w("python3 -I round2/20_delta_rerun.py \\")
    w("    --repo <repo root> \\")
    w("    --tested round2/tested_variants.csv \\")
    w("    --out round2/20_RERUN_DELTA.md")
    w("```")
    w("")
    w("It reads the pre-fix tables from git and the corrected ones from the")
    w("working tree, so it only works while the re-run is still uncommitted.")
    w("After the commit, compare against the parent: `git show <commit>^:<path>`.")
    w("")

    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote %s (%d tables)" % (a.out, len(reports)))


if __name__ == "__main__":
    main()

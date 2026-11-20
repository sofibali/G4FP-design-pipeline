#!/usr/bin/env python3
"""
13_delta_ipae_gate.py -- the SWITCH gate, rebuilt as a confidence contrast.

WHY THIS EXISTS
---------------
The round-2 SWITCH gate measures a chromophore-local apo/holo RMSD after
superposition.  round2/14_KNOWN_LIMITATIONS.md records that signal as
0.19-0.22 A at a superposition residual of 0.16-0.24 A: the quantity is
smaller than the fit error that produces it, so the gate ranks noise.

Chang & Polizzi (bioRxiv 2026, doi:10.64898/2026.08.03.742366) designed a
working small-molecule-driven conformational switch and did NOT gate on a
structural difference.  They gated on a *confidence contrast*: keep designs
where the structure predictor is confident about the mobile element WITH
ligand and uncertain about it WITHOUT.  Their hits had dPAE > 3.

This script computes the same quantity for G4FP, between two AlphaFold3
runs that already exist on disk.

    MODULE  = the redesigned G4-facing region  (chain A 15-60 + 240-285,
              i.e. `ligandmpnn.design_residues` from pipeline_config.yaml)
    CORE    = the rest of chain A, minus a sequence buffer either side of
              the module

    iPAE(state) = mean PAE over all MODULE<->CORE token pairs, both
                  directions (PAE is asymmetric)
    dPAE        = iPAE(apo) - iPAE(holo)

A design that switches is floppy at the interface without G4 and pinned with
it, so dPAE is positive and large.  A design that is rigid in both states --
which is what round 1 selected for -- scores ~0.

WHAT MAKES THIS DIFFERENT FROM THE METRIC IT REPLACES
-----------------------------------------------------
No superposition anywhere.  PAE is expressed in the model's own frame, so
the metric cannot encode its own fit error (defect 4 of the
scoring-metric-integrity checklist).  Nothing is re-predicted: both PAE
matrices are already in the `*_confidences.json` AlphaFold3 wrote.

Also emitted, from the same files and at no extra cost:
    bind_ipae_holo  mean PAE across the protein<->G4 interface (the BIND gate)
    bind_pae_min    chain_pair_pae_min[A][B] from the summary file
    dpae_contact    dPAE using a contact-derived module instead of the
                    design window, as an independent check that the result
                    does not depend on where the window was drawn

USAGE
    python3 13_delta_ipae_gate.py --root .. --out 13_delta_ipae.csv
    python3 13_delta_ipae_gate.py --root .. --limit 20     # smoke test
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import os
import statistics
import subprocess
import sys
from typing import Dict, List, Sequence, Tuple

# ----------------------------------------------------------------------------
# Definitions.  Keep these here, not in a config file, so the CSV header can
# record exactly what produced the numbers.
# ----------------------------------------------------------------------------

#: Redesigned, G4-facing region.  `ligandmpnn.design_residues` in
#: pipeline_config.yaml, identical across all ten templates (verified).
DESIGN_WINDOW: Tuple[Tuple[int, int], ...] = ((15, 60), (240, 285))

#: Residues within this many positions of a module boundary are dropped from
#: CORE.  PAE between sequence neighbours is low whatever the fold does, so
#: without a buffer the contrast is diluted by trivially-confident pairs.
SEQ_BUFFER = 4

#: Protein chain, G4 chain.  AlphaFold3 writes the protein first.
PROTEIN_CHAIN = "A"
LIGAND_CHAIN = "B"

#: Contact probability above which a protein token counts as G4-facing, for
#: the `dpae_contact` cross-check.  AF3 contact_probs, holo run only.
CONTACT_P = 0.5


#: Per-template MODULE definition, read from 13_fixed_core_definition.csv when
#: --module-map is given.  That file records, per design, where the grafted
#: RHAU DSM module actually sits (`rhau_N`, `rhau_C`) and where the barrel is.
#: The design window is NOT the module: on most templates the module is at
#: 1-28 + 279-306 while the window is 15-60 + 240-285, which is mostly barrel.
#: Scoring the window against "the rest" therefore measures barrel against
#: barrel.  Prefer --module-map; fall back to the window only to reproduce the
#: first run.
MODULE_MAP_CSV = "13_fixed_core_definition.csv"


def parse_ranges(spec: str) -> Tuple[Tuple[int, int], ...]:
    out = []
    for part in spec.replace(",", " ").split():
        lo, hi = part.split("-")
        out.append((int(lo), int(hi)))
    return tuple(out)


def load_module_map(path: str) -> Dict[str, Tuple[Tuple[Tuple[int, int], ...], Tuple[Tuple[int, int], ...]]]:
    """template -> (module ranges, core ranges), from the modal row per template.

    Templates whose rows disagree are returned anyway but reported, because a
    module whose position is not constant across designs of one template is a
    fact about that template, not a reason to drop it silently.
    """
    import collections

    seen = collections.defaultdict(collections.Counter)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            tpl = r["input"].split("_d0")[0]
            seen[tpl][(r["rhau_N"], r["rhau_C"], r["barrel"])] += 1
    out = {}
    for tpl, counter in seen.items():
        (n, c, barrel), _ = counter.most_common(1)[0]
        out[tpl] = (parse_ranges(n) + parse_ranges(c), parse_ranges(barrel))
        if len(counter) > 1:
            print(
                "  NOTE %-20s module position varies across its designs (%d distinct); using the modal one %s + %s"
                % (tpl, len(counter), n, c),
                file=sys.stderr,
            )
    return out


def in_ranges(res_id: int, ranges: Sequence[Tuple[int, int]], pad: int = 0) -> bool:
    return any(lo - pad <= res_id <= hi + pad for lo, hi in ranges)


def in_window(res_id: int) -> bool:
    return in_ranges(res_id, DESIGN_WINDOW)


def near_window(res_id: int) -> bool:
    return in_ranges(res_id, DESIGN_WINDOW, SEQ_BUFFER)


# ----------------------------------------------------------------------------
# Core computation
# ----------------------------------------------------------------------------


def mean_block(pae: Sequence[Sequence[float]], rows: Sequence[int], cols: Sequence[int]) -> float:
    """Mean PAE over rows x cols in both directions.

    PAE[i][j] is the expected error in residue j's position when the frame of
    residue i is aligned; it is not symmetric, so both blocks are averaged.
    """
    if not rows or not cols:
        raise ValueError("empty index set -- module/core definition selected nothing")
    total = 0.0
    n = 0
    for i in rows:
        row = pae[i]
        for j in cols:
            total += row[j]
            n += 1
    for i in cols:
        row = pae[i]
        for j in rows:
            total += row[j]
            n += 1
    return total / n


def load_conf(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def protein_index(conf: dict) -> Dict[int, int]:
    """residue id -> token index, for the protein chain only."""
    out: Dict[int, int] = {}
    for idx, (ch, res) in enumerate(zip(conf["token_chain_ids"], conf["token_res_ids"])):
        if ch != PROTEIN_CHAIN:
            continue
        if res in out:
            # A residue spanning several tokens would make a residue->token map
            # ambiguous.  For a polypeptide AF3 emits one token per residue; if
            # that ever stops being true, fail rather than silently pick one.
            raise ValueError("residue %d maps to more than one token" % res)
        out[res] = idx
    return out


def contact_module(conf: dict) -> List[int]:
    """Protein residues the holo model places in contact with the G4."""
    cp = conf["contact_probs"]
    chains = conf["token_chain_ids"]
    lig = [i for i, ch in enumerate(chains) if ch == LIGAND_CHAIN]
    if not lig:
        raise ValueError("no ligand chain %s in holo prediction" % LIGAND_CHAIN)
    out = []
    for res, idx in protein_index(conf).items():
        if max(cp[idx][j] for j in lig) >= CONTACT_P:
            out.append(res)
    return out


def score_design(holo_path: str, apo_path: str, mod_core=None) -> dict:
    holo = load_conf(holo_path)
    apo = load_conf(apo_path)

    h_idx = protein_index(holo)
    a_idx = protein_index(apo)

    # Fail loudly rather than scoring misaligned residues.  This is the check
    # that makes the comparison meaningful at all.
    if set(h_idx) != set(a_idx):
        raise ValueError(
            "apo/holo protein residue sets differ (holo %d, apo %d, symmetric "
            "difference %d) -- cannot compare states"
            % (len(h_idx), len(a_idx), len(set(h_idx) ^ set(a_idx)))
        )

    residues = sorted(h_idx)
    if mod_core is None:
        module = [r for r in residues if in_window(r)]
        core = [r for r in residues if not near_window(r)]
        label = "design window"
    else:
        mod_r, core_r = mod_core
        module = [r for r in residues if in_ranges(r, mod_r)]
        core = [r for r in residues if in_ranges(r, core_r) and not in_ranges(r, mod_r, SEQ_BUFFER)]
        label = "RHAU module vs barrel"
    if not module or not core:
        raise ValueError(
            "%s selects nothing on a %d-residue chain (module %d, core %d)"
            % (label, len(residues), len(module), len(core))
        )

    h_pae, a_pae = holo["pae"], apo["pae"]
    ipae_holo = mean_block(h_pae, [h_idx[r] for r in module], [h_idx[r] for r in core])
    ipae_apo = mean_block(a_pae, [a_idx[r] for r in module], [a_idx[r] for r in core])

    # Independent module definition: whatever the holo model actually touches.
    try:
        cmod = [r for r in contact_module(holo) if r in h_idx]
        ccore = [r for r in residues if r not in set(cmod)]
        # same sequence buffer, applied to the contact set
        cnear = set()
        for r in cmod:
            cnear.update(range(r - SEQ_BUFFER, r + SEQ_BUFFER + 1))
        ccore = [r for r in ccore if r not in cnear]
        if cmod and ccore:
            dpae_contact = mean_block(
                a_pae, [a_idx[r] for r in cmod], [a_idx[r] for r in ccore]
            ) - mean_block(h_pae, [h_idx[r] for r in cmod], [h_idx[r] for r in ccore])
            n_contact = len(cmod)
        else:
            dpae_contact, n_contact = float("nan"), len(cmod)
    except ValueError:
        dpae_contact, n_contact = float("nan"), 0

    # BIND: protein <-> G4 interface confidence, holo only.
    lig_tokens = [i for i, ch in enumerate(holo["token_chain_ids"]) if ch == LIGAND_CHAIN]
    bind_ipae = (
        mean_block(h_pae, [h_idx[r] for r in residues], lig_tokens)
        if lig_tokens
        else float("nan")
    )

    summ_path = holo_path.replace("_confidences.json", "_summary_confidences.json")
    bind_pae_min, iptm, ptm_holo = float("nan"), float("nan"), float("nan")
    if os.path.exists(summ_path):
        s = load_conf(summ_path)
        iptm = s.get("iptm") if s.get("iptm") is not None else float("nan")
        ptm_holo = s.get("ptm", float("nan"))
        cpm = s.get("chain_pair_pae_min")
        if cpm and len(cpm) > 1 and cpm[0][1] is not None:
            bind_pae_min = cpm[0][1]

    return {
        "n_residues": len(residues),
        "n_module_window": len(module),
        "n_core_window": len(core),
        "n_module_contact": n_contact,
        "ipae_apo": round(ipae_apo, 4),
        "ipae_holo": round(ipae_holo, 4),
        "dpae": round(ipae_apo - ipae_holo, 4),
        "dpae_contact": round(dpae_contact, 4) if dpae_contact == dpae_contact else "",
        "bind_ipae_holo": round(bind_ipae, 4) if bind_ipae == bind_ipae else "",
        "bind_pae_min": round(bind_pae_min, 4) if bind_pae_min == bind_pae_min else "",
        "iptm_holo": iptm if iptm == iptm else "",
        "ptm_holo": ptm_holo,
    }


# ----------------------------------------------------------------------------
# Walking the tree
# ----------------------------------------------------------------------------


def find_pairs(root: str) -> List[Tuple[str, str, str, str]]:
    """(template, design_id, holo_confidences, apo_confidences), sorted."""
    pairs = []
    for entry in sorted(os.listdir(root)):
        if not entry.startswith("output_"):
            continue
        bdir = os.path.join(root, entry, "03_alphafold3_predictions_bound")
        adir = os.path.join(root, entry, "03_alphafold3_predictions_apo")
        if not (os.path.isdir(bdir) and os.path.isdir(adir)):
            continue
        for sub in sorted(os.listdir(bdir)):
            if not sub.endswith("_bound"):
                continue
            did = sub[: -len("_bound")]
            holo = os.path.join(bdir, sub, "%s_bound_confidences.json" % did)
            apo = os.path.join(adir, "%s_apo" % did, "%s_apo_confidences.json" % did)
            if os.path.exists(holo) and os.path.exists(apo):
                pairs.append((entry, did, holo, apo))
    return pairs


def git_hash(root: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", root, "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unknown"


FIELDS = [
    "template",
    "design_id",
    "n_residues",
    "n_module_window",
    "n_core_window",
    "n_module_contact",
    "ipae_apo",
    "ipae_holo",
    "dpae",
    "dpae_contact",
    "bind_ipae_holo",
    "bind_pae_min",
    "iptm_holo",
    "ptm_holo",
]


def variance_report(rows: List[dict]) -> None:
    """Defect-3 check: a column with no variance contributes nothing.

    Run this every time.  A metric that is constant or mostly empty is not
    participating in any ranking built from it, however it is weighted.
    """
    print("\n-- column participation (non-empty / distinct / stdev) --", file=sys.stderr)
    for col in FIELDS[2:]:
        vals = []
        for r in rows:
            v = r[col]
            if v == "" or v is None:
                continue
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                pass
        if not vals:
            print("  %-18s 0 values -- NOT PARTICIPATING" % col, file=sys.stderr)
            continue
        sd = statistics.pstdev(vals) if len(vals) > 1 else 0.0
        flag = "  <-- zero variance, NOT PARTICIPATING" if sd == 0 else ""
        print(
            "  %-18s %5d / %5d distinct  sd=%8.4f%s"
            % (col, len(vals), len(set(vals)), sd, flag),
            file=sys.stderr,
        )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="..", help="campaign root holding output_* directories")
    ap.add_argument("--out", default="13_delta_ipae.csv")
    ap.add_argument("--limit", type=int, default=0, help="score only the first N pairs")
    ap.add_argument(
        "--module-map",
        nargs="?",
        const=MODULE_MAP_CSV,
        default=None,
        help="use the per-template RHAU module / barrel split from this CSV "
        "instead of the design window (recommended; default file: %s)" % MODULE_MAP_CSV,
    )
    args = ap.parse_args()

    mmap = None
    if args.module_map:
        path = args.module_map
        if not os.path.isabs(path):
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
        print("module definition: RHAU module vs barrel, from %s" % path, file=sys.stderr)
        mmap = load_module_map(path)
    else:
        print("module definition: design window %s (mostly barrel -- prefer --module-map)"
              % " + ".join("%d-%d" % w for w in DESIGN_WINDOW), file=sys.stderr)

    root = os.path.abspath(args.root)
    pairs = find_pairs(root)
    if args.limit:
        pairs = pairs[: args.limit]
    if not pairs:
        print("no apo/holo pairs found under %s" % root, file=sys.stderr)
        return 1
    print("scoring %d apo/holo pairs" % len(pairs), file=sys.stderr)

    rows: List[dict] = []
    failures: List[Tuple[str, str, str]] = []
    for n, (tpl, did, holo, apo) in enumerate(pairs, 1):
        try:
            key = tpl[len("output_G4FP_"):] if tpl.startswith("output_G4FP_") else tpl
            rec = score_design(holo, apo, mmap.get(key) if mmap else None)
            if mmap and key not in mmap:
                raise ValueError("no module definition for template %r" % key)
        except Exception as exc:  # recorded, never silently skipped
            failures.append((tpl, did, "%s: %s" % (type(exc).__name__, exc)))
            continue
        rec.update(template=tpl, design_id=did)
        rows.append(rec)
        if n % 100 == 0:
            print("  %d/%d" % (n, len(pairs)), file=sys.stderr)

    out = os.path.abspath(args.out)
    with open(out, "w", newline="") as fh:
        fh.write("# 13_delta_ipae_gate.py\n")
        fh.write("# generated: %s\n" % _dt.datetime.now().isoformat(timespec="seconds"))
        fh.write("# git: %s\n" % git_hash(root))
        fh.write("# root: %s\n" % root)
        if mmap:
            fh.write("# module (MODULE): grafted RHAU DSM segments, per template, from %s\n" % MODULE_MAP_CSV)
            fh.write("# core (CORE): the barrel, per template, from the same file\n")
            for k in sorted(mmap):
                fh.write("#   %-20s module %s | core %s\n"
                         % (k, " + ".join("%d-%d" % r for r in mmap[k][0]),
                            " + ".join("%d-%d" % r for r in mmap[k][1])))
        else:
            fh.write("# module (MODULE): chain %s residues %s  [pipeline_config.yaml ligandmpnn.design_residues]\n"
                     % (PROTEIN_CHAIN, " + ".join("%d-%d" % w for w in DESIGN_WINDOW)))
            fh.write("# core (CORE): remaining chain %s residues\n" % PROTEIN_CHAIN)
        fh.write("# both: core excludes +/-%d residues around each module boundary\n" % SEQ_BUFFER)
        fh.write("# dpae = mean PAE(MODULE<->CORE) apo - same holo; both directions; no superposition\n")
        fh.write("# contact module threshold: AF3 contact_probs >= %.2f against chain %s\n"
                 % (CONTACT_P, LIGAND_CHAIN))
        fh.write("# rows: %d   failures: %d\n" % (len(rows), len(failures)))
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})

    print("wrote %s  (%d rows, %d failures)" % (out, len(rows), len(failures)), file=sys.stderr)
    if failures:
        print("\n-- failures --", file=sys.stderr)
        for tpl, did, msg in failures[:20]:
            print("  %s/%s  %s" % (tpl, did, msg), file=sys.stderr)
        if len(failures) > 20:
            print("  ... and %d more" % (len(failures) - 20), file=sys.stderr)
    variance_report(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())

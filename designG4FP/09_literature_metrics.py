#!/usr/bin/env python3
"""Literature-supported design metrics for the round-1 G4FP designs.

Computed entirely from AF3 output already on disk -- no new prediction.

WHY THIS EXISTS
---------------
The round-1 fitness score was a hand-weighted sum of AlphaFold confidence
values.  Roughly half that weight sat on quantities the literature says AF
cannot report:

  * chromophore pLDDT (25%) -- chromophore formation is post-translational
    cyclization of the Thr-Tyr-Gly triad; AF emits only the pre-reaction
    backbone (Hartley et al., PLoS ONE 17:e0267560, 2022).
  * apo pLDDT vs template, INVERTED (10%) -- pLDDT difference does not track
    experimental ddG (Pak et al., PLoS ONE 18:e0282689, 2023), and AF3 is
    memorisation-biased toward the folded/holo state for over-represented
    folds (Yu et al., PNAS 123:e2530709123, 2026).
  * holo-vs-apo RMSD (15%) -- undercut by the same ligand-insensitivity.

This module replaces that composite with a small set of THRESHOLDED,
individually-validated measures, which is what the field actually uses
(Bennett et al., Nat Commun 14:2625, 2023; Watson et al., Nature 620:1089,
2023).  Nothing here is summed into a single score on purpose.

METRICS
-------
pae_interaction_ordered
    Mean interchain PAE between ORDERED protein residues (pLDDT >= 70) and
    the DNA.  This is the Bennett 2023 metric (`pAE_interaction < 10`), the
    single best-validated interface filter in de novo binder design.
    Restricting to ordered residues is required here: see below.

min_contact_distance / n_clash_lt_2_5 / clash_frac_terminal
    Steric sanity.  Every round-1 holo model contains protein-DNA heavy-atom
    contacts far closer than physically possible (minimum 0.30 A observed;
    a covalent C-O bond is 1.43 A).  94% of the residues involved lie in the
    first or last 30 positions and carry pLDDT ~55-69 against ~91 for the
    chain -- AF3 threading disordered termini through the duplex.  These
    contaminate any global interface score, ipTM included.

n_core_contacts
    Ordered-residue protein-DNA heavy-atom pairs in the real contact range
    (3.0-4.5 A).  A genuine interface has many; a spurious one has few.

chromophore_cyclization_distance
    Thr(N) to Gly(C) distance across the triad -- the bond that must form for
    the chromophore to mature.  A physical geometric criterion in the spirit
    of the catalytic-triad distance gates used in Choi et al., "De novo design
    of cysteine proteases", Nature (2026).  Reported INSTEAD of treating
    chromophore pLDDT as a fluorescence proxy.

chromophore_plddt / apo_mean_plddt / apo_fraction_disordered
    Reported as descriptors, NOT inverted and NOT scored.  AF3 predicted a
    confidently folded apo state for every round-1 design; these columns
    record that fact rather than pretending to rank on it.
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from chromophore_utils import triad_start_in_sequence, THREE2ONE

ORDERED_PLDDT = 70.0      # Bennett-style: ignore residues AF itself calls unreliable
CLASH_CUTOFF = 2.5        # A, heavy-atom; below this is not a contact
CONTACT_LO, CONTACT_HI = 3.0, 4.5   # A, real non-bonded contact range
TERMINAL_WINDOW = 30      # residues at each end


def load_cif(cif: Path):
    """Return (cols, rows) of the _atom_site loop."""
    cols: Dict[str, int] = {}
    rows: List[List[str]] = []
    in_loop = False
    with open(cif) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('_atom_site.'):
                cols[s.split('.', 1)[1]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if s.startswith(('ATOM', 'HETATM')):
                    rows.append(s.split())
                elif s.startswith('#') and rows:
                    break
    return cols, rows


def _col(cols, *names):
    for n in names:
        if n in cols:
            return cols[n]
    raise KeyError(names)


def structure_metrics(cif: Path, want_dna: bool) -> Optional[Dict]:
    cols, rows = load_cif(cif)
    if not rows:
        return None
    cx, cy, cz = cols['Cartn_x'], cols['Cartn_y'], cols['Cartn_z']
    ca = _col(cols, 'auth_asym_id', 'label_asym_id')
    cs = _col(cols, 'auth_seq_id', 'label_seq_id')
    cc = _col(cols, 'label_comp_id', 'auth_comp_id')
    cat = _col(cols, 'label_atom_id', 'auth_atom_id')
    cb = cols['B_iso_or_equiv']
    el = cols['type_symbol']

    prot = [f for f in rows if f[ca] == 'A' and f[el] != 'H']
    if not prot:
        return None
    dna = [f for f in rows if f[ca] == 'B' and f[el] != 'H']

    # per-residue pLDDT (B-factor column) and sequence
    acc = defaultdict(list)
    comp: Dict[int, str] = {}
    for f in prot:
        sid = int(f[cs])
        acc[sid].append(float(f[cb]))
        comp.setdefault(sid, f[cc])
    order = sorted(acc)
    plddt = {s: float(np.mean(acc[s])) for s in order}
    seq = ''.join(THREE2ONE.get(comp[s], 'X') for s in order)
    nres = len(order)

    out: Dict[str, float] = {
        'n_residues': nres,
        'mean_plddt': float(np.mean(list(plddt.values()))),
        'frac_residues_ordered': float(np.mean([plddt[s] >= ORDERED_PLDDT for s in order])),
    }

    # ---- chromophore: corrected index + physical cyclization geometry -----
    t0 = triad_start_in_sequence(seq)          # 0-based position in `order`
    if t0 is not None and t0 + 2 < nres:
        triad = [order[t0], order[t0 + 1], order[t0 + 2]]
        out['chromophore_start'] = t0 + 1
        out['chromophore_plddt'] = float(np.mean([plddt[s] for s in triad]))
        # Thr N  ->  Gly C : the bond formed on cyclization
        def xyz(sid, atom):
            for f in prot:
                if int(f[cs]) == sid and f[cat] == atom:
                    return np.array([float(f[cx]), float(f[cy]), float(f[cz])])
            return None
        n_thr, c_gly = xyz(triad[0], 'N'), xyz(triad[2], 'C')
        if n_thr is not None and c_gly is not None:
            out['chromophore_cyclization_distance'] = float(np.linalg.norm(n_thr - c_gly))

    if not want_dna or not dna:
        return out

    # ---- interface: steric sanity and ordered-residue contacts ------------
    P = np.array([[float(f[cx]), float(f[cy]), float(f[cz])] for f in prot])
    D = np.array([[float(f[cx]), float(f[cy]), float(f[cz])] for f in dna])
    pres = np.array([int(f[cs]) for f in prot])
    d = np.linalg.norm(P[:, None, :] - D[None, :, :], axis=2)
    dmin_per_atom = d.min(axis=1)

    out['min_contact_distance'] = float(d.min())
    clash_atoms = np.where(dmin_per_atom < CLASH_CUTOFF)[0]
    clash_res = sorted({int(pres[i]) for i in clash_atoms})
    out['n_clash_lt_2_5'] = int(len(clash_atoms))
    out['n_clash_residues'] = int(len(clash_res))
    if clash_res:
        term = sum(1 for r in clash_res if r <= TERMINAL_WINDOW or r > nres - TERMINAL_WINDOW)
        out['clash_frac_terminal'] = term / len(clash_res)
        out['clash_median_plddt'] = float(np.median([plddt[r] for r in clash_res]))

    ordered_mask = np.array([plddt.get(r, 0.0) >= ORDERED_PLDDT for r in pres])
    if ordered_mask.any():
        d_ord = d[ordered_mask]
        in_range = (d_ord >= CONTACT_LO) & (d_ord <= CONTACT_HI)
        out['n_core_contacts'] = int(in_range.sum())
        out['min_contact_distance_ordered'] = float(d_ord.min())
        out['n_clash_lt_2_5_ordered'] = int((d_ord < CLASH_CUTOFF).sum())
    return out


def interface_pae_ordered(conf_json: Path, plddt_by_pos: np.ndarray) -> Optional[float]:
    """Mean interchain PAE between ordered protein tokens and the DNA tokens."""
    try:
        with open(conf_json) as fh:
            fc = json.load(fh)
    except Exception:
        return None
    if 'pae' not in fc or 'token_chain_ids' not in fc:
        return None
    pae = np.array(fc['pae'])
    tch = np.array(fc['token_chain_ids'])
    if pae.ndim != 2:
        return None
    prot_idx = np.where(tch == 'A')[0]
    other_idx = np.where(tch != 'A')[0]
    if len(prot_idx) == 0 or len(other_idx) == 0:
        return None
    n = min(len(prot_idx), len(plddt_by_pos))
    keep = prot_idx[:n][plddt_by_pos[:n] >= ORDERED_PLDDT]
    if len(keep) == 0:
        return None
    return float(np.mean(pae[np.ix_(keep, other_idx)]))


def per_residue_plddt_array(cif: Path) -> Optional[np.ndarray]:
    cols, rows = load_cif(cif)
    if not rows:
        return None
    ca = _col(cols, 'auth_asym_id', 'label_asym_id')
    cs = _col(cols, 'auth_seq_id', 'label_seq_id')
    cb, el = cols['B_iso_or_equiv'], cols['type_symbol']
    acc = defaultdict(list)
    for f in rows:
        if f[ca] == 'A' and f[el] != 'H':
            acc[int(f[cs])].append(float(f[cb]))
    if not acc:
        return None
    return np.array([float(np.mean(acc[s])) for s in sorted(acc)])


SEED_RE = re.compile(r'^seed-(\d+)_sample-(\d+)$')


def seed_dirs(design_dir: Path, limit: int) -> List[Path]:
    out = [d for d in sorted(design_dir.iterdir()) if d.is_dir() and SEED_RE.match(d.name)]
    return out[:limit]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--max-seeds', type=int, default=5,
                    help='seed/sample dirs per design per state (default 5)')
    ap.add_argument('--output', '-o', default='09_literature_metrics.csv')
    ap.add_argument('--templates', default='output_G4FP_*')
    args = ap.parse_args()

    records = []
    for od in sorted(glob.glob(args.templates)):
        tpl = Path(od).name.replace('output_', '')
        hroot = Path(od) / '03_alphafold3_predictions_bound'
        aroot = Path(od) / '03_alphafold3_predictions_apo'
        if not hroot.exists():
            continue
        designs = sorted(d for d in hroot.iterdir()
                         if d.is_dir() and re.match(r'^design_\d+_bound$', d.name))
        print(f"{tpl}: {len(designs)} designs")
        for hd in designs:
            m = re.match(r'^design_(\d+)_bound$', hd.name)
            seq_id = int(m.group(1))
            ad = aroot / f'design_{seq_id}_apo'
            rec = {'template': tpl, 'seq_id': seq_id,
                   'global_id': f'{tpl}_design_{seq_id:04d}'}

            for state, root in (('holo', hd), ('apo', ad)):
                if not root.exists():
                    continue
                vals = defaultdict(list)
                for sd in seed_dirs(root, args.max_seeds):
                    cif = sd / 'model.cif'
                    if not cif.exists():
                        continue
                    sm = structure_metrics(cif, want_dna=(state == 'holo'))
                    if sm:
                        for k, v in sm.items():
                            vals[k].append(v)
                    if state == 'holo':
                        pl = per_residue_plddt_array(cif)
                        if pl is not None:
                            ip = interface_pae_ordered(sd / 'confidences.json', pl)
                            if ip is not None:
                                vals['pae_interaction_ordered'].append(ip)
                for k, v in vals.items():
                    rec[f'{state}_{k}'] = float(np.mean(v))
                    if len(v) > 1:
                        rec[f'{state}_{k}_sd'] = float(np.std(v, ddof=1))
                rec[f'{state}_n_seeds'] = len(seed_dirs(root, args.max_seeds)) if root.exists() else 0
            records.append(rec)

    df = pd.DataFrame(records)
    df.to_csv(args.output, index=False)
    print(f"\nwrote {args.output}: {len(df)} designs x {len(df.columns)} columns")


if __name__ == '__main__':
    main()

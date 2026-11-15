"""Shared chromophore location and per-residue pLDDT extraction.

Both 05_analyze_af3_structures.py and 06_compare_ligand_states.py previously
carried their own copies of these two functions, and both copies were wrong in
the same two ways.  They live here now so there is one implementation to fix.

Defect 1 -- per-residue pLDDT by uniform atom binning
    The old code estimated residue boundaries as len(atoms) / n_residues and
    sliced the atom array into equal bins.  Real residues hold 4-14 atoms, so
    the mapping drifted cumulatively.  On a 306-residue chain the slice
    labelled "residues 197-199" read atoms belonging to residues 200-203,
    reporting pLDDT 90.6 where the true chromophore value was 73.9.

Defect 2 -- chromophore located by PDB residue NUMBER
    The old detect_chromophore_range() read the CRO heteroresidue's residue
    number out of the template PDB and used it as an index into AF3's
    sequential per-residue array.  CRO is numbered 197 in all ten templates,
    so the "auto-detect" returned (197, 199) every time -- identical to the
    hardcode it replaced.  But the templates carry numbering gaps before the
    chromophore, so its SEQUENTIAL position is 187-197 depending on the
    template, and AF3 output is renumbered sequentially.  Six of ten templates
    had their chromophore pLDDT measured on unrelated barrel residues.

Together these fed holo_chromophore_plddt (25% of the round-1 fitness weight)
and chromophore_plddt_diff (5%), inflated by +13 to +20 pLDDT units.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}

# The chromophore-forming tripeptide sits in an invariant cpGFP context.
# Anchoring on that context is the only reliable way to locate it, because
# residue NUMBERING differs between the template PDBs (numbering gaps) and
# AF3 output (renumbered sequentially).  A bare "TYG" search is not enough:
# a spurious TYG occurs ~170 residues into every one of these designs.
CHROMOPHORE_MOTIF = re.compile(r"VTTL(TYG)VQ")
CHROMOPHORE_MOTIF_LOOSE = re.compile(r"TTL(TYG)V")

_UNIFORM_BIN_WARNED = set()


def parse_cif_atoms(cif_path: Path) -> List[Tuple[str, int, str, str, float]]:
    """Parse an AF3 _model.cif into (chain, seq_id, comp_id, atom_id, bfactor).

    AF3 writes the per-atom pLDDT into the B-factor column.  Column positions
    are read from the _atom_site loop header rather than assumed, so this
    survives column reordering between AF3 versions.
    """
    atoms: List[Tuple[str, int, str, str, float]] = []
    cols: Dict[str, int] = {}
    in_loop = False
    with open(cif_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('_atom_site.'):
                cols[s.split('.', 1)[1]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if s.startswith(('ATOM', 'HETATM')):
                    f = s.split()
                    try:
                        chain = f[cols.get('auth_asym_id', cols.get('label_asym_id'))]
                        sid = int(f[cols.get('auth_seq_id', cols.get('label_seq_id'))])
                        comp = f[cols.get('auth_comp_id', cols.get('label_comp_id'))]
                        aid = f[cols.get('auth_atom_id', cols.get('label_atom_id'))]
                        b = float(f[cols['B_iso_or_equiv']])
                    except (KeyError, IndexError, ValueError, TypeError):
                        continue
                    atoms.append((chain, sid, comp, aid, b))
                elif s.startswith('#') and atoms:
                    break
    return atoms


def chain_sequence(atoms, chain_id: str = 'A') -> Tuple[str, List[int]]:
    """One-letter sequence of a chain, plus its residue ids in order."""
    comp: Dict[int, str] = {}
    for ch, sid, cp, _aid, _b in atoms:
        if ch == chain_id:
            comp.setdefault(sid, cp)
    order = sorted(comp)
    seq = ''.join(THREE2ONE.get(comp[s], 'X') for s in order)
    return seq, order


def find_model_cif(d: Path) -> Optional[Path]:
    """Any AF3 model .cif at or below *d* (first match).

    AF3 names these inconsistently: seed/sample directories hold a plain
    ``model.cif`` while top-level design directories hold
    ``<jobname>_model.cif``.  The glob must cover both -- matching only
    ``*_model.cif`` silently misses every seed directory and sends this
    module down its unreliable fallback path.
    """
    for pattern in ("model.cif", "*model.cif"):
        hits = sorted(d.glob(pattern))
        if hits:
            return hits[0]
    hits = sorted(d.rglob("*model.cif"))
    return hits[0] if hits else None


def triad_start_in_sequence(seq: str) -> Optional[int]:
    """0-based index of the T of the TYG triad in *seq*, or None.

    A bare "TYG" search is NOT safe: every one of these designs carries a
    spurious TYG roughly 170 residues in.  The surrounding cpGFP context is
    what makes the match reliable.
    """
    for pattern in (CHROMOPHORE_MOTIF, CHROMOPHORE_MOTIF_LOOSE):
        m = pattern.search(seq)
        if m:
            return m.start(1)
    return None


def detect_chromophore_range(output_dir: Path) -> Tuple[int, int]:
    """Locate the Thr-Tyr-Gly triad in AF3 SEQUENTIAL numbering.

    Returns 1-indexed sequential positions (start, end) for use as
    ``per_residue_array[start-1 : end]``.

    Raises ValueError rather than falling back to a hardcoded range: a silent
    fallback is what produced the round-1 error.
    """
    cif = find_model_cif(output_dir)
    if cif is not None:
        seq, _order = chain_sequence(parse_cif_atoms(cif), 'A')
        for pattern in (CHROMOPHORE_MOTIF, CHROMOPHORE_MOTIF_LOOSE):
            m = pattern.search(seq)
            if m:
                cs = m.start(1) + 1
                print(f"  Chromophore: TYG at sequential {cs}-{cs + 2} "
                      f"(motif {pattern.pattern}, from {cif.name})")
                return (cs, cs + 2)

    # Fallback: CRO's ORDINAL position in the template chain, not its number.
    tpl_name = output_dir.name.replace("output_", "")
    tpl_pdb = output_dir.parent / "inputs" / f"{tpl_name}.pdb"
    if tpl_pdb.exists():
        seen: List[str] = []
        cro_ordinal = None
        with open(tpl_pdb) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")) and line[21] == 'A':
                    key = line[22:27]
                    if not seen or seen[-1] != key:
                        seen.append(key)
                        if line[17:20].strip() == "CRO":
                            cro_ordinal = len(seen)
        if cro_ordinal is not None:
            print(f"  Chromophore: CRO at ordinal {cro_ordinal} in "
                  f"{tpl_pdb.name} (no AF3 model found; motif preferred)")
            return (cro_ordinal, cro_ordinal + 2)

    raise ValueError(
        f"Could not locate the chromophore triad for {output_dir.name}. "
        f"Refusing to fall back to a hardcoded range -- that silent fallback "
        f"is what put 30% of the round-1 fitness weight on the wrong residues."
    )


def compute_per_residue_plddt(full_conf: Dict, chain_id: str = 'A',
                              model_dir: Optional[Path] = None
                              ) -> Optional[np.ndarray]:
    """Mean pLDDT per residue for *chain_id*, in sequential residue order.

    Grouped by ACTUAL residue id from the model .cif.  Pass *model_dir* to get
    correct values; without it this falls back to the old uniform-bin estimate
    and warns, because those per-residue values are not reliable.
    """
    if model_dir is not None:
        cif = find_model_cif(Path(model_dir))
        if cif is not None:
            try:
                acc: Dict[int, List[float]] = defaultdict(list)
                for ch, sid, _comp, _aid, b in parse_cif_atoms(cif):
                    if ch == chain_id and not np.isnan(b):
                        acc[sid].append(b)
                if acc:
                    return np.array([float(np.mean(acc[s])) for s in sorted(acc)])
            except Exception:
                pass

    if 'atom_plddts' not in full_conf or 'atom_chain_ids' not in full_conf:
        return None
    atom_plddts = np.array(full_conf['atom_plddts'])
    atom_chains = full_conf['atom_chain_ids']
    chain_mask = np.array([ch == chain_id for ch in atom_chains])
    chain_plddts = atom_plddts[chain_mask]
    if len(chain_plddts) == 0:
        return None
    if 'token_chain_ids' in full_conf:
        token_chains = full_conf['token_chain_ids']
        n_residues = sum(1 for ch in token_chains if ch == chain_id)
        if n_residues > 0:
            key = str(model_dir)
            if key not in _UNIFORM_BIN_WARNED:
                _UNIFORM_BIN_WARNED.add(key)
                print(f"  WARNING: no model .cif for {model_dir}; falling back "
                      f"to uniform atom binning. Per-residue values -- "
                      f"chromophore pLDDT above all -- are NOT reliable.")
            atoms_per_residue = len(chain_plddts) / n_residues
            per_res = []
            for i in range(n_residues):
                start = int(round(i * atoms_per_residue))
                end = int(round((i + 1) * atoms_per_residue))
                if start < len(chain_plddts):
                    per_res.append(np.mean(chain_plddts[start:end]))
            return np.array(per_res)
    return chain_plddts

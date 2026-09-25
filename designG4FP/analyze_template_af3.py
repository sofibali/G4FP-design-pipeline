#!/usr/bin/env python3
"""
Analyze template AF3 predictions (25 seeds x bound + apo) and save
metrics to a CSV that 06 and 07 can read.

This avoids re-running the full 06 pipeline just for templates.

Input:  template_af3_outputs/{bound,apo}/template_*/ with seed dirs
Output: template_af3_metrics.csv

Usage:
    python3 analyze_template_af3.py
"""

import json
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List


def discover_seed_dirs(result_dir: Path) -> List[Path]:
    """Find seed-N_sample-M directories."""
    pattern = re.compile(r'^seed-\d+_sample-\d+$')
    return sorted(d for d in result_dir.iterdir()
                  if d.is_dir() and pattern.match(d.name))


def find_result_dir(base_dir: Path) -> Optional[Path]:
    """Find the dir containing seed subdirectories."""
    if not base_dir.exists():
        return None
    # Check for seeds directly
    if any(base_dir.glob("seed-*_sample-*")):
        return base_dir
    # Check one level down
    for sub in sorted(base_dir.iterdir()):
        if sub.is_dir() and any(sub.glob("seed-*_sample-*")):
            return sub
    return None


def analyze_seed(seed_dir: Path, chain_id: str = 'A') -> Optional[Dict]:
    """Extract metrics from one seed-sample directory."""
    result = {}

    # Summary confidences
    for sf in [seed_dir / 'summary_confidences.json',
               *seed_dir.glob('*_summary_confidences.json')]:
        if sf.exists():
            with open(sf) as f:
                sc = json.load(f)
            for k in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered']:
                if k in sc:
                    result[k] = sc[k]
            break

    # Full confidences (pLDDT, PAE)
    conf_files = [f for f in [seed_dir / 'confidences.json',
                               *seed_dir.glob('*_confidences.json')]
                  if f.exists() and 'summary' not in f.name]
    if conf_files:
        with open(conf_files[0]) as f:
            fc = json.load(f)

        # Per-residue pLDDT for chain
        plddts = fc.get('atom_plddts', [])
        chain_ids = fc.get('atom_chain_ids', [])
        if plddts and chain_ids:
            chain_plddts = [p for p, c in zip(plddts, chain_ids) if c == chain_id]
            if chain_plddts:
                result['mean_plddt'] = float(np.mean(chain_plddts))

                # Chromophore pLDDT (residues 197-199, TYG)
                # Estimate atoms per residue to extract chromophore
                token_chains = fc.get('token_chain_ids', [])
                if token_chains:
                    n_residues = sum(1 for ch in token_chains if ch == chain_id)
                    if n_residues > 0:
                        atoms_per_res = len(chain_plddts) / n_residues
                        # Chromophore at residues 197-199 (0-indexed: 196-198)
                        cs, ce = 196, 199  # 0-indexed start, exclusive end
                        if n_residues >= ce:
                            chrom_atoms = []
                            for ri in range(cs, ce):
                                start = int(round(ri * atoms_per_res))
                                end = int(round((ri + 1) * atoms_per_res))
                                chrom_atoms.extend(chain_plddts[start:end])
                            if chrom_atoms:
                                result['chromophore_plddt'] = float(np.mean(chrom_atoms))

        # PAE
        pae = fc.get('pae', [])
        if pae:
            pae_arr = np.array(pae)
            result['mean_pae'] = float(np.mean(pae_arr))

            # Interface PAE (protein vs DNA)
            if token_chains:
                n_prot = sum(1 for ch in token_chains if ch == chain_id)
                if pae_arr.ndim == 2 and pae_arr.shape[0] > n_prot:
                    result['interface_pae'] = float(
                        np.mean(pae_arr[:n_prot, n_prot:]))

    return result if result else None


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=str, default=None,
                        help='Path to template AF3 outputs (default: template_af3_outputs/)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output CSV path (default: template_af3_metrics.csv)')
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent

    if args.input_dir:
        tpl_af3_dir = Path(args.input_dir)
    else:
        tpl_af3_dir = base_dir / "template_af3_outputs"

    if not tpl_af3_dir.exists():
        print(f"ERROR: {tpl_af3_dir} not found")
        sys.exit(1)

    records = []

    for disk_state in ['bound', 'apo']:
        display_state = 'holo' if disk_state == 'bound' else 'apo'
        state_dir = tpl_af3_dir / disk_state
        if not state_dir.exists():
            continue

        for tpl_dir in sorted(state_dir.iterdir()):
            if not tpl_dir.is_dir():
                continue

            result_dir = find_result_dir(tpl_dir)
            if result_dir is None:
                continue

            seed_dirs = discover_seed_dirs(result_dir)
            if not seed_dirs:
                continue

            # Extract template name
            tname = tpl_dir.name
            for suffix in [f'_{disk_state}']:
                if tname.lower().endswith(suffix):
                    tname = tname[:-len(suffix)]
            if tname.lower().startswith('template_'):
                tname = tname[len('template_'):]

            # Normalize to match design template names (G4FP_...)
            # Handle lowercase dirs: g4fp_des1_cro_mod0 -> G4FP_des1_cro_mod0
            if tname.startswith('g4fp_') or tname.startswith('G4FP_'):
                tname = 'G4FP_' + tname[5:]

            # Analyze all seeds
            metrics = []
            for sd in seed_dirs:
                m = analyze_seed(sd)
                if m:
                    metrics.append(m)

            if not metrics:
                continue

            # Find or create record for this template
            existing = None
            for r in records:
                if r['template_name'] == tname:
                    existing = r
                    break
            if existing is None:
                existing = {'template_name': tname}
                records.append(existing)

            # Aggregate across seeds
            df_m = pd.DataFrame(metrics)
            prefix = f'{display_state}_'
            for col in df_m.columns:
                vals = df_m[col].dropna()
                if len(vals) > 0:
                    existing[f'{prefix}{col}'] = float(vals.mean())
                    existing[f'{prefix}{col}_sd'] = (
                        float(vals.std()) if len(vals) > 1 else 0.0)

            existing[f'n_{display_state}_seeds'] = len(metrics)

            print(f"  {tname} / {display_state}: {len(metrics)} seeds analyzed")

    if not records:
        print("ERROR: No template AF3 results found")
        sys.exit(1)

    # ── Compute matched-seed holo-vs-apo RMSD ──
    print("\nComputing holo-vs-apo RMSD for templates...")
    try:
        from Bio.PDB import MMCIFParser, Superimposer
        cif_parser = MMCIFParser(QUIET=True)

        def _extract_ca(cif_path):
            struct = cif_parser.get_structure('s', str(cif_path))
            ca = []
            for model in struct:
                for chain in model:
                    if chain.get_id() == 'A':
                        for res in chain:
                            if res.has_id('CA'):
                                ca.append(res['CA'])
            return ca

        for rec in records:
            tname = rec['template_name']
            # Find holo and apo seed dirs
            holo_seeds = apo_seeds = None
            for disk_state in ['bound', 'apo']:
                state_dir = tpl_af3_dir / disk_state
                if not state_dir.exists():
                    continue
                for tpl_dir in state_dir.iterdir():
                    if not tpl_dir.is_dir():
                        continue
                    rd = find_result_dir(tpl_dir)
                    if rd is None:
                        continue
                    sds = discover_seed_dirs(rd)
                    if not sds:
                        continue
                    # Check name matches
                    dn = tpl_dir.name
                    for sfx in [f'_{disk_state}']:
                        if dn.lower().endswith(sfx):
                            dn = dn[:-len(sfx)]
                    if dn.lower().startswith('template_'):
                        dn = dn[len('template_'):]
                    if dn.startswith('g4fp_') or dn.startswith('G4FP_'):
                        dn = 'G4FP_' + dn[5:]
                    if dn != tname:
                        continue
                    if disk_state == 'bound':
                        holo_seeds = sds
                    else:
                        apo_seeds = sds
                    break

            if not holo_seeds or not apo_seeds:
                continue

            # Match seeds by index (min of both lists)
            n_match = min(len(holo_seeds), len(apo_seeds))
            rmsds = []
            chrom_rmsds = []
            cs, ce = 196, 199  # 0-indexed chromophore range

            for i in range(n_match):
                h_cif = holo_seeds[i] / 'model.cif'
                a_cif = apo_seeds[i] / 'model.cif'
                if not h_cif.exists() or not a_cif.exists():
                    continue
                try:
                    h_ca = _extract_ca(h_cif)
                    a_ca = _extract_ca(a_cif)
                    min_len = min(len(h_ca), len(a_ca))
                    if min_len < 10:
                        continue
                    h_ca = h_ca[:min_len]
                    a_ca = a_ca[:min_len]
                    sup = Superimposer()
                    sup.set_atoms(h_ca, a_ca)
                    rmsds.append(sup.rms)
                    if min_len >= ce:
                        sup_c = Superimposer()
                        sup_c.set_atoms(h_ca[cs:ce], a_ca[cs:ce])
                        chrom_rmsds.append(sup_c.rms)
                except Exception:
                    continue

            if rmsds:
                rec['global_rmsd'] = float(np.mean(rmsds))
                rec['global_rmsd_sd'] = float(np.std(rmsds)) if len(rmsds) > 1 else 0.0
            if chrom_rmsds:
                rec['chromophore_rmsd'] = float(np.mean(chrom_rmsds))
                rec['chromophore_rmsd_sd'] = float(np.std(chrom_rmsds)) if len(chrom_rmsds) > 1 else 0.0

            print(f"  {tname}: global_rmsd={rec.get('global_rmsd', 'N/A'):.2f}A "
                  f"({n_match} matched seed pairs)")

    except ImportError:
        print("  WARNING: BioPython not available, skipping RMSD")

    # ── Compute difference metrics ──
    for rec in records:
        hp = rec.get('holo_mean_plddt')
        ap = rec.get('apo_mean_plddt')
        if hp is not None and ap is not None:
            rec['mean_plddt_diff'] = hp - ap
            hsd = rec.get('holo_mean_plddt_sd', 0)
            asd = rec.get('apo_mean_plddt_sd', 0)
            rec['mean_plddt_diff_sd'] = float(np.sqrt(hsd**2 + asd**2))

        hcp = rec.get('holo_chromophore_plddt')
        acp = rec.get('apo_chromophore_plddt')
        if hcp is not None and acp is not None:
            rec['chromophore_plddt_diff'] = hcp - acp

    # ── Compute fitness score (same weights as 07) ──
    df = pd.DataFrame(records)

    def _norm(series):
        smin, smax = series.min(), series.max()
        if smax == smin:
            return pd.Series(0.5, index=series.index)
        return (series - smin) / (smax - smin)

    # Components
    if 'holo_chromophore_plddt' in df.columns:
        df['_n_chrom_holo'] = _norm(df['holo_chromophore_plddt'].fillna(0))
    else:
        df['_n_chrom_holo'] = 0.5

    if 'holo_iptm' in df.columns:
        df['_n_iptm'] = _norm(df['holo_iptm'].fillna(0))
    else:
        df['_n_iptm'] = 0.5

    if 'global_rmsd' in df.columns:
        df['_n_rmsd'] = _norm(df['global_rmsd'].fillna(0))
    else:
        df['_n_rmsd'] = 0.0

    if 'holo_ptm' in df.columns:
        df['_n_holo_ptm'] = _norm(df['holo_ptm'].fillna(0))
    else:
        df['_n_holo_ptm'] = 0.0

    # For templates, rel_apo = 0 (they are the reference)
    df['_n_rel_apo'] = 0.5

    if 'mean_plddt_diff_sd' in df.columns:
        df['_n_confidence'] = 1.0 - _norm(df['mean_plddt_diff_sd'].fillna(0))
    else:
        df['_n_confidence'] = 0.5

    if 'chromophore_plddt_diff' in df.columns:
        df['_n_chrom_diff'] = _norm(df['chromophore_plddt_diff'].fillna(0))
    else:
        df['_n_chrom_diff'] = 0.0

    df['fitness_score'] = (
        0.25 * df['_n_chrom_holo']
        + 0.25 * df['_n_iptm']
        + 0.15 * df['_n_rmsd']
        + 0.15 * df['_n_holo_ptm']
        + 0.10 * df['_n_rel_apo']
        + 0.05 * df['_n_confidence']
        + 0.05 * df['_n_chrom_diff']
    )

    # Drop internal columns
    df = df.drop(columns=[c for c in df.columns if c.startswith('_n_')])

    # Sort columns for readability
    id_cols = ['template_name']
    score_cols = ['fitness_score']
    other_cols = sorted([c for c in df.columns if c not in id_cols + score_cols])
    df = df[id_cols + score_cols + other_cols]

    if args.output:
        out_path = Path(args.output)
    else:
        out_path = base_dir / "template_af3_metrics.csv"
    df.to_csv(out_path, index=False)

    print(f"\nSaved: {out_path}")
    print(f"  {len(df)} templates, {len(df.columns)} columns")

    # Print summary
    print(f"\nTemplate summary:")
    for _, row in df.iterrows():
        tpl = row['template_name']
        fs = row.get('fitness_score', float('nan'))
        hp = row.get('holo_mean_plddt', float('nan'))
        ap = row.get('apo_mean_plddt', float('nan'))
        hi = row.get('holo_iptm', float('nan'))
        rmsd = row.get('global_rmsd', float('nan'))
        print(f"  {tpl}: fitness={fs:.3f}  holo_plddt={hp:.1f}  "
              f"apo_plddt={ap:.1f}  iPTM={hi:.3f}  RMSD={rmsd:.2f}A")


if __name__ == '__main__':
    main()

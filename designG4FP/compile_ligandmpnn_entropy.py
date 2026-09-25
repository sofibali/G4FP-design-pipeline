#!/usr/bin/env python3
"""
Compile LigandMPNN entropy metrics across all templates into CSV files.

Outputs:
  07_results/07_ligandmpnn_entropy.csv           — per-design entropy + confidence
  07_results/07_ligandmpnn_positional_entropy.csv — aligned grid: template sequences
                                                    as rows, alignment columns,
                                                    positional entropy per template

Usage:
    python3 compile_ligandmpnn_entropy.py
"""

import numpy as np
import pandas as pd
import torch
from pathlib import Path
from itertools import combinations
import sys


AA = 'ACDEFGHIKLMNPQRSTVWY'


def compute_entropy(probs, eps=1e-10):
    """Shannon entropy from probability distribution (base 2, bits)."""
    p = np.clip(probs, eps, 1.0)
    return -np.sum(p * np.log2(p), axis=-1)


def decode_sequence(indices):
    """Convert integer indices to amino acid string."""
    return ''.join(AA[i] if i < 20 else 'X' for i in indices)


def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Simple NW alignment. Returns aligned sequences with gaps."""
    n, m = len(seq1), len(seq2)
    # Score matrix
    dp = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        dp[i][0] = dp[i-1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j-1] + gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1] + s, dp[i-1][j] + gap, dp[i][j-1] + gap)

    # Traceback
    aln1, aln2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + s:
                aln1.append(seq1[i-1])
                aln2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            aln1.append(seq1[i-1])
            aln2.append('-')
            i -= 1
        else:
            aln1.append('-')
            aln2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(aln1)), ''.join(reversed(aln2))


def progressive_msa(sequences, names):
    """Simple progressive MSA using pairwise NW alignment.

    Aligns all sequences to the longest one, then merges into a single alignment.
    """
    # Use longest sequence as reference
    ref_idx = max(range(len(sequences)), key=lambda i: len(sequences[i]))
    ref_seq = sequences[ref_idx]

    # Align each sequence to reference
    aligned = {}
    ref_aligned_versions = {}
    for i, (seq, name) in enumerate(zip(sequences, names)):
        if i == ref_idx:
            aligned[name] = seq
            ref_aligned_versions[name] = seq
            continue
        aln_ref, aln_seq = needleman_wunsch(ref_seq, seq)
        aligned[name] = aln_seq
        ref_aligned_versions[name] = aln_ref

    # Merge: the reference alignment may have different gap patterns
    # when aligned to different sequences. Build a master alignment
    # by inserting all gaps.

    # Start with the reference as the master
    master_ref = list(ref_seq)
    master_positions = list(range(len(ref_seq)))  # original ref positions

    # For each other sequence, find where it inserts gaps in the reference
    all_insertions = {}  # position -> set of names that need insertions
    for name in names:
        if name == names[ref_idx]:
            continue
        aln_ref = ref_aligned_versions[name]
        ref_pos = -1
        insertions = []
        for col_idx, (rc, sc) in enumerate(zip(aln_ref, aligned[name])):
            if rc != '-':
                ref_pos += 1
            else:
                # Gap in reference = insertion in this sequence
                insertions.append((ref_pos, col_idx))
        all_insertions[name] = insertions

    # Build master alignment: for simplicity, align each to ref independently
    # and pad with gaps to the max alignment length
    pairwise_alns = {}
    pairwise_ref_alns = {}
    for name in names:
        if name == names[ref_idx]:
            pairwise_alns[name] = ref_seq
            pairwise_ref_alns[name] = ref_seq
        else:
            aln_ref, aln_seq = needleman_wunsch(ref_seq, sequences[names.index(name)])
            pairwise_alns[name] = aln_seq
            pairwise_ref_alns[name] = aln_ref

    # Map each alignment back to reference coordinates
    # Build column mapping: for each pairwise alignment, map columns to ref positions
    result_columns = []  # list of dicts: {name: (residue_char, original_position)}

    # Use the longest pairwise alignment as the column template
    max_aln_len = max(len(pairwise_alns[n]) for n in names)

    # Simpler approach: use reference as backbone, insert columns for each unique gap
    # For each pairwise alignment with ref, track ref_col -> aln mapping
    ref_to_aln = {}
    for name in names:
        aln_ref = pairwise_ref_alns[name]
        aln_seq = pairwise_alns[name]
        mapping = []
        ref_col = -1
        for aln_col, (rc, sc) in enumerate(zip(aln_ref, aln_seq)):
            if rc != '-':
                ref_col += 1
                mapping.append(('match', ref_col, sc))
            else:
                mapping.append(('insert', ref_col, sc))
        ref_to_aln[name] = mapping

    # Build final alignment: iterate through ref positions,
    # collecting insertions between each pair of ref positions
    final_aln = {name: [] for name in names}
    final_pos_map = {name: [] for name in names}  # original position index or -1

    for ref_pos in range(len(ref_seq)):
        # Add the match column
        for name in names:
            found = False
            for entry in ref_to_aln[name]:
                if entry[0] == 'match' and entry[1] == ref_pos:
                    final_aln[name].append(entry[2])
                    # Find original position
                    orig_pos = sum(1 for e in ref_to_aln[name][:ref_to_aln[name].index(entry)+1]
                                   if e[2] != '-') - 1
                    final_pos_map[name].append(orig_pos if entry[2] != '-' else -1)
                    found = True
                    break
            if not found:
                final_aln[name].append('-')
                final_pos_map[name].append(-1)

    # Convert to strings
    result = {}
    pos_maps = {}
    for name in names:
        result[name] = ''.join(final_aln[name])
        pos_maps[name] = final_pos_map[name]

    return result, pos_maps


def main():
    base_dir = Path(__file__).resolve().parent

    # ── Load all template data ──
    templates = {}
    for odir in sorted(base_dir.glob("output_G4FP_*")):
        tpl_name = odir.name.replace("output_", "")
        stats_dir = odir / "01_ligandmpnn" / "stats"
        if not stats_dir.exists():
            continue
        pt_files = list(stats_dir.glob("*.pt"))
        if not pt_files:
            continue

        print(f"  Loading {tpl_name}...")
        data = torch.load(pt_files[0], map_location='cpu', weights_only=False)

        sampling_probs = data['sampling_probs'].numpy()  # (N, L, 20)
        native_seq = data['native_sequence'].numpy()
        chain_mask = data['chain_mask'].numpy()

        seq_str = decode_sequence(native_seq)
        entropy_per_pos = compute_entropy(sampling_probs)  # (N, L)
        mean_entropy = np.mean(entropy_per_pos, axis=0)    # (L,)
        designed_mask = chain_mask > 0.5

        templates[tpl_name] = {
            'sequence': seq_str,
            'mean_entropy': mean_entropy,
            'designed_mask': designed_mask,
            'entropy_matrix': entropy_per_pos,
            'sampling_probs': sampling_probs,
            'chain_mask': chain_mask,
        }

    if not templates:
        print("ERROR: No stats files found")
        sys.exit(1)

    names = sorted(templates.keys())
    sequences = [templates[n]['sequence'] for n in names]

    # ── Per-design entropy CSV (same as before) ──
    print("\nCompiling per-design entropy...")
    all_records = []
    for tpl_name in names:
        tdata = templates[tpl_name]
        ent_matrix = tdata['entropy_matrix']
        designed = tdata['designed_mask']
        n_designs = ent_matrix.shape[0]

        # Chromophore region (±5 around 197-199, 0-indexed)
        n_res = len(tdata['sequence'])
        chrom_region = np.zeros(n_res, dtype=bool)
        cs = max(0, 191)  # 196-5, 0-indexed
        ce = min(n_res, 204)  # 199+5, 0-indexed
        chrom_region[cs:ce] = True

        # Load filtered sequences for confidence scores
        filt_csv = base_dir / f"output_{tpl_name}" / "02_filtered_sequences.csv"
        conf_map = {}
        id_to_seq = {}
        if filt_csv.exists():
            filt_df = pd.read_csv(filt_csv)
            if 'id' in filt_df.columns and 'seq_id' in filt_df.columns:
                for _, row in filt_df.iterrows():
                    sid = int(row['seq_id'])
                    did = int(row['id'])
                    id_to_seq[did] = sid
                    conf_map[sid] = {
                        'overall_confidence': row.get('overall_confidence'),
                        'ligand_confidence': row.get('ligand_confidence'),
                        'seq_rec': row.get('seq_rec'),
                        'combined_confidence': row.get('combined_confidence'),
                    }

        for i in range(n_designs):
            if i not in id_to_seq:
                continue
            sid = id_to_seq[i]
            ent = ent_matrix[i]
            rec = {
                'global_id': f'{tpl_name}_design_{sid:04d}',
                'template': tpl_name,
                'seq_id': sid,
                'mean_entropy_all': float(np.mean(ent)),
                'mean_entropy_designed': float(np.mean(ent[designed])) if designed.any() else np.nan,
                'mean_entropy_fixed': float(np.mean(ent[~designed])) if (~designed).any() else np.nan,
                'mean_entropy_chromophore_region': float(np.mean(ent[chrom_region])) if chrom_region.any() else np.nan,
                'max_entropy': float(np.max(ent)),
                'min_entropy_designed': float(np.min(ent[designed])) if designed.any() else np.nan,
                'n_designed_positions': int(designed.sum()),
                'n_total_positions': n_res,
            }
            if sid in conf_map:
                rec.update(conf_map[sid])
            all_records.append(rec)

        print(f"    {tpl_name}: {sum(1 for r in all_records if r['template'] == tpl_name)} designs")

    df_designs = pd.DataFrame(all_records)
    results_dir = base_dir / "07_results"
    results_dir.mkdir(exist_ok=True, parents=True)
    df_designs.to_csv(results_dir / "07_ligandmpnn_entropy.csv", index=False)
    print(f"  Saved 07_ligandmpnn_entropy.csv ({len(df_designs)} designs)")

    # ── Positional entropy grid CSV ──
    print("\nBuilding aligned positional entropy grid...")

    # Align all sequences to longest
    ref_idx = max(range(len(sequences)), key=lambda i: len(sequences[i]))
    ref_seq = sequences[ref_idx]
    ref_name = names[ref_idx]

    # Pairwise align each to reference
    alignments = {}
    orig_pos_maps = {}  # name -> list of original 0-indexed positions (-1 for gap)

    for i, (seq, name) in enumerate(zip(sequences, names)):
        if i == ref_idx:
            alignments[name] = seq
            orig_pos_maps[name] = list(range(len(seq)))
            continue
        aln_ref, aln_seq = needleman_wunsch(ref_seq, seq)

        # Now we need to reconcile: the reference alignment from this pair
        # may have insertions. We'll handle this by noting that the reference
        # is the longest and these are circular permutants of the same protein,
        # so insertions in the reference should be minimal.

        # Map alignment back to reference coordinates
        ref_pos = -1
        seq_pos = -1
        mapped_seq = []
        mapped_orig = []

        for rc, sc in zip(aln_ref, aln_seq):
            if rc != '-':
                ref_pos += 1
            if sc != '-':
                seq_pos += 1

            if rc != '-':
                # This column maps to a reference position
                mapped_seq.append(sc)
                mapped_orig.append(seq_pos if sc != '-' else -1)
            # If rc == '-', this is an insertion in the query relative to ref
            # We skip it (these are rare for same-protein circular permutants)

        alignments[name] = ''.join(mapped_seq)
        orig_pos_maps[name] = mapped_orig

    # All alignments should now be the same length (= ref length)
    aln_len = len(ref_seq)

    # Build the grid
    # Row 1: alignment column number (1-indexed)
    # Row 2: reference sequence
    # Rows 3+: each template with sequence and entropy

    grid_data = {'alignment_column': list(range(1, aln_len + 1))}

    for name in names:
        aln_seq = alignments[name]
        orig_map = orig_pos_maps[name]
        mean_ent = templates[name]['mean_entropy']
        designed = templates[name]['designed_mask']

        # Pad if shorter than ref
        while len(aln_seq) < aln_len:
            aln_seq += '-'
            orig_map.append(-1)

        seq_col = []
        ent_col = []
        designed_col = []

        for col in range(aln_len):
            seq_col.append(aln_seq[col])
            orig_idx = orig_map[col] if col < len(orig_map) else -1
            if orig_idx >= 0 and orig_idx < len(mean_ent):
                ent_col.append(round(float(mean_ent[orig_idx]), 4))
                designed_col.append(bool(designed[orig_idx]))
            else:
                ent_col.append(np.nan)
                designed_col.append(False)

        short = name.replace('G4FP_', '')
        grid_data[f'{short}_residue'] = seq_col
        grid_data[f'{short}_entropy'] = ent_col
        grid_data[f'{short}_designed'] = designed_col

    df_grid = pd.DataFrame(grid_data)

    pos_path = results_dir / "07_ligandmpnn_positional_entropy.csv"
    df_grid.to_csv(pos_path, index=False)
    print(f"  Saved 07_ligandmpnn_positional_entropy.csv "
          f"({aln_len} columns x {len(names)} templates)")

    # Summary
    print(f"\nEntropy summary by template:")
    for name in names:
        designed = templates[name]['designed_mask']
        ent = templates[name]['mean_entropy']
        des_ent = ent[designed]
        print(f"  {name}: designed_entropy={des_ent.mean():.3f} "
              f"(min={des_ent.min():.3f}, max={des_ent.max():.3f})")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Step 8: Export final G4FP candidates for gene synthesis ordering.

Reads 07_results/07_final_candidates.csv and produces formats suitable for:
  - Gene synthesis ordering (IDT, Twist, GenScript)
  - 96-well plate mapping
  - Lab notebook summary

Outputs in 08_synthesis_order/:
  08_synthesis_sequences.csv    -- sequences with names for ordering
  08_plate_map.csv              -- 96-well plate layout
  08_lab_summary.csv            -- compact summary for lab notebook
  08_sequences.fa               -- clean FASTA for synthesis vendor upload

Usage:
    python 08_export_for_synthesis.py
    python 08_export_for_synthesis.py --n-plates 3 --prefix G4FP_v1
    python 08_export_for_synthesis.py --add-flanks --5prime ATGGTG --3prime TGATAA
"""

import sys
import argparse
import string
from pathlib import Path

import pandas as pd

# Standard codon table (E. coli optimized, most common codons)
CODON_TABLE_ECOLI = {
    "A": "GCG", "R": "CGT", "N": "AAC", "D": "GAT", "C": "TGC",
    "E": "GAA", "Q": "CAG", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCG",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAT", "V": "GTG",
    "*": "TAA",
}

# Human-optimized codon table
CODON_TABLE_HUMAN = {
    "A": "GCC", "R": "AGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "E": "GAG", "Q": "CAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


def reverse_translate(aa_seq: str, table: dict) -> str:
    """Convert amino acid sequence to DNA using a codon table."""
    return "".join(table.get(aa, "NNN") for aa in aa_seq.upper())


def generate_plate_map(n_sequences: int, n_plates: int = 0):
    """
    Generate 96-well plate positions.

    Returns list of (plate, well) tuples like ("Plate1", "A01").
    """
    rows = list(string.ascii_uppercase[:8])  # A-H
    cols = list(range(1, 13))                # 1-12
    wells_per_plate = 96

    if n_plates <= 0:
        n_plates = (n_sequences + wells_per_plate - 1) // wells_per_plate

    positions = []
    for plate_idx in range(n_plates):
        for row in rows:
            for col in cols:
                if len(positions) >= n_sequences:
                    return positions
                positions.append((f"Plate{plate_idx + 1}", f"{row}{col:02d}"))

    return positions


def main():
    parser = argparse.ArgumentParser(
        description="Export G4FP candidates for gene synthesis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=str, default=None,
                        help="Input CSV (default: 07_results/07_final_candidates.csv)")
    parser.add_argument("--prefix", type=str, default="G4FP",
                        help="Construct name prefix (default: G4FP)")
    parser.add_argument("--organism", type=str, default="ecoli",
                        choices=["ecoli", "human"],
                        help="Codon optimization target (default: ecoli)")
    parser.add_argument("--add-flanks", action="store_true",
                        help="Add 5' and 3' flanking sequences")
    parser.add_argument("--5prime", type=str, default="",
                        dest="flank_5",
                        help="5' flanking DNA sequence (e.g., start codon region)")
    parser.add_argument("--3prime", type=str, default="",
                        dest="flank_3",
                        help="3' flanking DNA sequence (e.g., stop codon, tag)")
    parser.add_argument("--add-stop", action="store_true", default=True,
                        help="Append stop codon (default: True)")
    parser.add_argument("--no-stop", action="store_false", dest="add_stop",
                        help="Do not append stop codon")
    parser.add_argument("--n-plates", type=int, default=0,
                        help="Number of 96-well plates (default: auto)")
    parser.add_argument("--max-sequences", type=int, default=0,
                        help="Limit number of sequences to export (0=all)")
    parser.add_argument("--no-dna", action="store_true",
                        help="Skip reverse translation, export protein only")

    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent
    output_dir = base_dir / "08_synthesis_order"
    output_dir.mkdir(exist_ok=True, parents=True)

    # Load candidates
    if args.input:
        input_csv = Path(args.input)
    else:
        input_csv = base_dir / "07_results" / "07_final_candidates.csv"

    if not input_csv.exists():
        print(f"ERROR: Input file not found: {input_csv}")
        print("       Run 07_aggregate_and_rank.py first.")
        sys.exit(1)

    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} candidates from {input_csv}")

    if args.max_sequences > 0:
        df = df.head(args.max_sequences)
        print(f"  Limited to top {len(df)}")

    if "sequence" not in df.columns:
        print("ERROR: No 'sequence' column found in input CSV.")
        sys.exit(1)

    # Select codon table
    codon_table = CODON_TABLE_ECOLI if args.organism == "ecoli" else CODON_TABLE_HUMAN

    # Generate plate positions
    positions = generate_plate_map(len(df), args.n_plates)

    # Build export table
    rows = []
    for i, (_, row) in enumerate(df.iterrows()):
        plate, well = positions[i] if i < len(positions) else ("overflow", f"X{i:03d}")

        construct_name = f"{args.prefix}_{i + 1:04d}"
        aa_seq = row["sequence"]
        aa_len = len(aa_seq)

        # Reverse translate
        if not args.no_dna:
            dna_seq = reverse_translate(aa_seq, codon_table)
            if args.add_stop:
                dna_seq += codon_table.get("*", "TAA")
            if args.add_flanks:
                dna_seq = args.flank_5 + dna_seq + args.flank_3
            dna_len = len(dna_seq)
        else:
            dna_seq = ""
            dna_len = 0

        entry = {
            "construct_name": construct_name,
            "plate": plate,
            "well": well,
            "protein_sequence": aa_seq,
            "protein_length": aa_len,
        }

        if not args.no_dna:
            entry["dna_sequence"] = dna_seq
            entry["dna_length"] = dna_len

        # Carry over key metrics
        for col in ["global_id", "template", "final_rank", "fitness_score",
                     "mean_plddt_diff", "bound_mean_plddt", "apo_mean_plddt",
                     "bound_ptm", "chromophore_plddt_diff"]:
            if col in row:
                entry[col] = row[col]

        rows.append(entry)

    export_df = pd.DataFrame(rows)

    # Save synthesis order CSV
    export_df.to_csv(output_dir / "08_synthesis_sequences.csv", index=False)
    print(f"\n  08_synthesis_sequences.csv: {len(export_df)} constructs")

    # Save plate map (compact)
    plate_cols = ["construct_name", "plate", "well", "protein_length"]
    if "fitness_score" in export_df.columns:
        plate_cols.append("fitness_score")
    export_df[plate_cols].to_csv(output_dir / "08_plate_map.csv", index=False)
    print(f"  08_plate_map.csv")

    # Lab summary
    summary_cols = ["construct_name", "plate", "well", "protein_length",
                    "template", "final_rank", "fitness_score"]
    summary_cols = [c for c in summary_cols if c in export_df.columns]
    for metric in ["mean_plddt_diff", "bound_mean_plddt", "apo_mean_plddt",
                    "bound_ptm", "chromophore_plddt_diff"]:
        if metric in export_df.columns:
            summary_cols.append(metric)
    export_df[summary_cols].to_csv(output_dir / "08_lab_summary.csv", index=False)
    print(f"  08_lab_summary.csv")

    # Clean FASTA
    with open(output_dir / "08_sequences.fa", "w") as f:
        for _, row in export_df.iterrows():
            f.write(f">{row['construct_name']}\n{row['protein_sequence']}\n")
    print(f"  08_sequences.fa")

    # DNA FASTA (if generated)
    if not args.no_dna:
        with open(output_dir / "08_dna_sequences.fa", "w") as f:
            for _, row in export_df.iterrows():
                f.write(f">{row['construct_name']}\n{row['dna_sequence']}\n")
        print(f"  08_dna_sequences.fa")

    # Summary stats
    print(f"\n{'=' * 60}")
    print(f"Synthesis Order Summary")
    print(f"{'=' * 60}")
    print(f"  Total constructs:  {len(export_df)}")
    print(f"  Plates needed:     {export_df['plate'].nunique()}")
    print(f"  Protein lengths:   {export_df['protein_length'].min()}-{export_df['protein_length'].max()} aa")
    if not args.no_dna:
        print(f"  DNA lengths:       {export_df['dna_length'].min()}-{export_df['dna_length'].max()} bp")
        print(f"  Codon table:       {args.organism}")
    if "template" in export_df.columns:
        print(f"  Templates:         {export_df['template'].nunique()}")
    print(f"\n  Output: {output_dir}/")
    print(f"\n  NOTE: For production orders, consider using a dedicated codon")
    print(f"  optimization tool (IDT Codon Opt, GenSmart, etc.) on the protein")
    print(f"  sequences. The reverse translation here uses most-common codons")
    print(f"  but does not optimize for GC content, repeats, or restriction sites.")


if __name__ == "__main__":
    main()

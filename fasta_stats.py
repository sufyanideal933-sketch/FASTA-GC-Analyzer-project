#!/usr/bin/env python3
"""
FASTA sequence statistics and GC-content plots.

Calculates sequence length, A/C/G/T counts, GC percentage, ambiguous-base
counts, and gap counts. Writes a CSV report and two PNG plots.

GC% = (G + C) / (A + C + G + T) * 100.
Ambiguous symbols (for example N, R, Y) and gaps are excluded from the
GC denominator. Sequence length includes all accepted symbols, including gaps.

Usage:
    python fasta_stats.py example.fasta
    python fasta_stats.py my_sequences.fasta --outdir results
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Standard IUPAC nucleotide ambiguity codes, plus common alignment gap symbols.
VALID_DNA = set("ACGTRYSWKMBDHVN-.")
BASES = "ACGT"


def read_fasta(path):
    """Parse FASTA records into a list of (sequence_id, sequence) pairs.

    Keeps record order and detects duplicate IDs rather than silently
    overwriting an earlier sequence. The ID is the first word after '>'.
    """
    records = []
    current_id = None
    chunks = []

    def save_current():
        if current_id is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence:
            raise ValueError(f"FASTA record {current_id!r} has an empty sequence.")
        invalid = sorted(set(sequence) - VALID_DNA)
        if invalid:
            chars = ", ".join(repr(char) for char in invalid)
            raise ValueError(
                f"FASTA record {current_id!r} contains unsupported character(s): {chars}. "
                "Accepted symbols are IUPAC DNA codes plus '-' and '.'."
            )
        records.append((current_id, sequence))

    with open(path, "r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                save_current()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header on line {line_number}.")
                current_id = header.split()[0]
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(
                        f"Sequence data found before the first FASTA header "
                        f"(line {line_number}). The file must start with a '>' header."
                    )
                # Remove internal spaces/tabs that may occur in formatted sequence lines.
                chunks.append("".join(line.split()))

    save_current()

    if not records:
        raise ValueError("No FASTA sequences found in the input file.")

    ids = [seq_id for seq_id, _ in records]
    duplicates = sorted({seq_id for seq_id in ids if ids.count(seq_id) > 1})
    if duplicates:
        raise ValueError("Duplicate FASTA sequence ID(s): " + ", ".join(duplicates))

    return records


def gc_content(seq):
    """Return (GC percentage, number of A/C/G/T bases used in denominator)."""
    denominator = sum(seq.count(base) for base in BASES)
    if denominator == 0:
        return float("nan"), 0
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / denominator, denominator


def build_table(records):
    """Build a DataFrame with per-sequence length, nucleotide counts, and GC%."""
    rows = []
    for seq_id, seq in records:
        gc_percent, gc_denominator = gc_content(seq)
        base_counts = {f"{base}_count": seq.count(base) for base in BASES}
        gap_count = seq.count("-") + seq.count(".")
        ambiguous_count = sum(seq.count(base) for base in set(seq) if base not in BASES + "-.")
        rows.append({
            "id": seq_id,
            "length": len(seq),
            **base_counts,
            "gc_percent": round(gc_percent, 2) if math.isfinite(gc_percent) else np.nan,
            "gc_denominator_acgt": gc_denominator,
            "ambiguous_count": ambiguous_count,
            "gap_count": gap_count,
            "n_count": seq.count("N"),
        })
    return pd.DataFrame(rows)


def make_plots(df, outdir):
    """Save GC histogram and length-vs-GC scatter plot, ignoring undefined GC%."""
    valid_gc = df.dropna(subset=["gc_percent"])
    if valid_gc.empty:
        print("Warning: no sequence contains A/C/G/T bases; GC plots were not created.",
              file=sys.stderr)
        return

    mean_gc = valid_gc["gc_percent"].mean()

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(valid_gc["gc_percent"], bins=min(10, max(1, len(valid_gc))),
            edgecolor="white")
    ax.axvline(mean_gc, linestyle="--", label=f"Mean = {mean_gc:.2f}%")
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Number of sequences")
    ax.set_title("GC content distribution")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "gc_histogram.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(valid_gc["length"], valid_gc["gc_percent"], alpha=0.8)
    ax.set_xlabel("Sequence length (symbols, including gaps)")
    ax.set_ylabel("GC content (%)")
    ax.set_title("Sequence length vs GC content")
    fig.tight_layout()
    fig.savefig(outdir / "length_vs_gc.png", dpi=200)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate FASTA sequence lengths and GC content, export CSV, and plot results."
    )
    parser.add_argument("fasta", help="path to an input FASTA file")
    parser.add_argument("--outdir", default="results", help="output folder (default: results)")
    args = parser.parse_args(argv)

    fasta_path = Path(args.fasta)
    if not fasta_path.is_file():
        parser.error(f"input FASTA file does not exist or is not a file: {fasta_path}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        records = read_fasta(fasta_path)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    df = build_table(records)
    csv_path = outdir / "sequence_stats.csv"
    df.to_csv(csv_path, index=False)
    make_plots(df, outdir)

    lengths = df["length"].to_numpy(dtype=float)
    valid_gc = df["gc_percent"].dropna()
    print(f"Sequences read     : {len(df)}")
    print(f"Length (symbols)   : min {int(lengths.min())}, max {int(lengths.max())}, "
          f"mean {np.mean(lengths):.1f}")
    print(f"Sequences with GC% : {len(valid_gc)} / {len(df)}")
    if len(valid_gc):
        std_gc = valid_gc.std(ddof=1) if len(valid_gc) > 1 else 0.0
        print(f"GC content (%)     : mean {valid_gc.mean():.2f}, std {std_gc:.2f}")
    print(f"CSV report         : {csv_path}")
    if len(valid_gc):
        print(f"GC histogram       : {outdir / 'gc_histogram.png'}")
        print(f"Length-vs-GC plot  : {outdir / 'length_vs_gc.png'}")


if __name__ == "__main__":
    main()

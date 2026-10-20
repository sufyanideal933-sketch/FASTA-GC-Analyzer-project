# FASTA GC Analyzer

A small command-line bioinformatics tool that reads a DNA FASTA file and reports per-sequence statistics.

## Features

- Parses multi-line FASTA records and preserves record order.
- Reports sequence length, A/C/G/T counts, GC percentage, GC denominator, ambiguous-base count, gap count, and N count.
- Validates sequence symbols and reports empty sequences, malformed headers, and duplicate IDs rather than silently producing misleading output.
- Exports a CSV table and creates a GC-content histogram and length-vs-GC scatter plot.
- Handles sequences with no A/C/G/T bases: GC% is reported as undefined (`NaN`) and such records are omitted from GC plots.

## GC definition

GC percentage is calculated as:

`(G + C) / (A + C + G + T) * 100`

Ambiguous IUPAC symbols (such as `N`, `R`, or `Y`) and alignment gaps (`-` or `.`) are excluded from the denominator. Sequence length includes all accepted symbols, including gaps. The accepted sequence alphabet is the standard IUPAC DNA alphabet plus `-` and `.`.

## Installation

Requires Python 3.8 or later.

```bash
pip install -r requirements.txt
```

## Usage

```bash
python fasta_stats.py example.fasta
python fasta_stats.py my_sequences.fasta --outdir my_results
```

Output files (in `results/` by default):

| File | Contents |
|---|---|
| `sequence_stats.csv` | One row per FASTA record with length, base counts, GC%, ambiguous-base count, gap count, and N count |
| `gc_histogram.png` | Distribution of defined GC percentages |
| `length_vs_gc.png` | Sequence length against GC percentage |

## Tests

Run the built-in unit tests from this directory:

```bash
python -m unittest -v
```

The tests cover multi-line/lowercase FASTA input, GC calculation with ambiguous symbols, invalid characters, duplicate IDs, empty records, and per-base counts.

## Example data

`example.fasta` contains synthetic, randomly generated sequences for software testing. They are not real biological sequences. For real analyses, use a FASTA file from a trusted source such as NCBI or Ensembl and retain its accession/source information.

## Interpretation notes

- This tool computes simple sequence composition statistics; it does not identify genes, variants, or biological function.
- A sequence containing only ambiguous bases or gaps has no defined GC percentage.
- Length includes gap symbols; when comparing aligned sequences, interpret lengths accordingly.

## License

MIT

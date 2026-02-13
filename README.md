# motif-mark

**motif-mark** is a visualization tool for exploring sequence motifs involved in eukaryotic alternative splicing.
Given a FASTA file and a list of sequence motifs, the script identifies where each motif occurs within exons and introns and produces a visual diagram highlighting their positions.

This tool is useful for quickly examining regulatory elements such as exonic splicing enhancers (ESEs), intronic splicing silencers (ISSs), and other short binding motifs that influence exon inclusion.

## Features

- Identifies motifs in DNA sequences provided as FASTA.
- Supports up to ten user‑defined motifs.
- Highlights motifs on exons and introns using distinct colors.
- Generates clear, annotated visualizations of motif distribution.

## Input Requirements

- FASTA File
    - May contain one or multiple sequences.
    - Exons should be indicated using uppercase; introns using lowercase.
        - Example: ATGAGTactgtgAACGT

- Motif File
    - A simple text file with one motif per line.
    - Motifs are case‑insensitive.
    - IUPAC degenerate bases are supported (Y, N, -, X, etc.)

## Using motif-mark

python motif-mark-oop.py -f <your_fasta_file> -m <motif_file>
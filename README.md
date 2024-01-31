# CirculoCov: Circular-Aware Coverage for Draft Genomes

## Overview
CirculoCov is a Python tool designed for circular-aware coverage analysis of draft genomes. This tool is designed to take a draft genome, map reads to it, and then visualize the coverage across the genome. This tool can also include the coverage of Illumina reads (paired-end or single-end).

## Requirements
- Python 3.9
- pandas
- matplotlib

## Installation
```
pip install circulocov
```

## Usage
```
circulocov -g fasta -n nanopore.fastq -i illumina1.fastq illumina2.fastq -o out
```

Options
-n, --nanopore: Path to the Nanopore fastq file (required).
-i, --illumina: Path to the Illumina fastq files (optional).
-o, --out: directory for results
-w, --window: window size for coverage

## Output
The output is
- A csv file with each contig broken into windows with their corresponding depths for Illumina and nanopore files
- png files showing depth

## Examples


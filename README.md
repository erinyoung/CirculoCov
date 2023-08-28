# WARNING : THIS TOOL IS UNDER DEVELOPMENT AND IS NOT ACTUALLY FUNCTIONAL

# Draft Genome Coverage

Draft Genome Covarege (DCG) is a tool to estimate and visualize the genome coverage of nanopore and/or illumina reads on a draft assembly.

Dependencies:
- pygenome-viz
- pip3

(Expected) usage:

```
dgc -n nanopore_reads.fastq -i illumina_R1.fastq illumina_R2.fastq -r draft_genome.fasta -o output_prefix
```

Input fastq and fasta files can be '.gz' compressed.

This will create a subset of files with the prefix set in the command line.

Output files
- dgc_coverage.png is a visualization of the coverage for each contig (not recommended for draft genomes with >10 contigs)
- dgc_coverage.csv is a csv file with a row for each contig and a column for contig name, percent coverage for nanopore reads, average depth for nanopore reads, percent coverage for illumina reads, average coverage for illumina reads

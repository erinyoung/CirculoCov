#!/usr/bin/env python3
# pylint: disable=logging-fstring-interpolation
# pylint: disable=R0912
# pylint: disable=R0915

'''
Author: Erin Young

Description:

This script will take
- reference
- reads
- map them with minimap2
- determine coverage
- visualize

EXAMPLE:
circulocov -g fasta -i illumina.fastq -n nanopore.fastq -o out
'''

import argparse
import logging
import os
import sys
import tempfile
import pandas as pd

#from utils.circular import circular
from circulocov.utils.genome_prep       import genome_prep
from circulocov.utils.mapping           import mapping
from circulocov.utils.counts            import counts
from circulocov.utils.merge_dataframe   import merge_cov_dataframe, merge_depth_dataframe
from circulocov.utils.summary           import summary
from circulocov.utils.visualize         import visualize
from circulocov.utils.extract           import extract

def main():
    """ Get coverage for draft genomes """

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 0. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    version = '0.1.20230104'

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample',
                        required = False,
                        type = str,
                        help = 'Sample name',
                        default= "circulocov")
    parser.add_argument('-g', '--genome',
                        required = True,
                        type = str,
                        help = 'Genome (draft or complete)')
    parser.add_argument('-i', '--illumina',
                        nargs = '+',
                        required = False,
                        type = str,
                        help = 'Input Illumina fastq(s)')
    parser.add_argument('-n', '--nanopore',
                        required = False,
                        type = str,
                        help ='Input nanopore fastq')
    parser.add_argument('-p', '--pacbio',
                        required = False,
                        type = str,
                        help ='Input pacbio fastq')
    parser.add_argument('-a', '--all',
                        action=argparse.BooleanOptionalAction)
    parser.add_argument('-d', '--padding',
                        required = False,
                        type = int,
                        help = 'Amount of padding added to circular sequences',
                        default = 10000)
    parser.add_argument('-w', '--window',
                        required = False,
                        type = int,
                        help = 'Window size for coverage',
                        default = 500)
    parser.add_argument('-o', '--out',
                        required = False,
                        type = str,
                        help = 'Result directory',
                        default = 'CirculoCov')
    parser.add_argument('-log', '--loglevel',
                        required = False,
                        type = str,
                        help = 'Logging level',
                        default = 'INFO')
    parser.add_argument('-t', '--threads',
                        required = False,
                        type = int,
                        help = 'Number of threads to use',
                        default=4)
    parser.add_argument('-v', '--version',
                        help='Print version and exit',
                        action = 'version',
                        version = version)
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s - %(message)s',
        datefmt = '%y-%b-%d %H:%M:%S',
        level=args.loglevel.upper())

    logging.info(f"Filenames :\t{args.sample}")
    logging.info(f"CirculoCov ver :\t{str(version)}")
    logging.info(f"Final directory :\t{str(args.out)}")
    logging.info(f"Genome file :\t{str(args.genome)}")
    logging.info(f"Num threads :\t{str(args.threads)}")
    logging.info(f"Padding length :\t{str(args.padding)}")
    if args.all:
        logging.info("All is set :\tWill create windows and graph coverages")
        logging.info(f"Window size :\t{str(args.window)}")
    if args.nanopore:
        logging.info(f"Nanopore file :\t{str(args.nanopore)}")
    if args.illumina:
        logging.info(f"Illumina file(s) :\t{', '.join(args.illumina)}")
    if args.pacbio:
        logging.info(f"PacBio file :\t{str(args.pacbio)}")

    if not args.nanopore and not args.illumina and not args.pacbio:
        logging.fatal('Cannot run without fastq files!')
        sys.exit(1)

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    with tempfile.TemporaryDirectory(dir = args.out) as temp_dir:

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 1. Setup                 #####
        ##### ----- ----- ----- ----- ----- #####

        genome_dict, fasta, total_length = genome_prep(args, temp_dir)

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 2. Minimap2              #####
        ##### ----- ----- ----- ----- ----- #####

        if args.nanopore:
            nanopore_bam = mapping(args.nanopore, fasta, "map-ont", args, temp_dir)
        else:
            nanopore_bam = ""

        if args.illumina:
            illumina_bam = mapping(args.illumina, fasta, "sr", args, temp_dir)
        else:
            illumina_bam = ""

        if args.pacbio:
            pacbio_bam   = mapping(args.pacbio, fasta, "map-pb", args, temp_dir)
        else:
            pacbio_bam   = ""

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 3. Samtools Counts       #####
        ##### ----- ----- ----- ----- ----- #####

        # used for testing
        # nanopore_bam = "CirculoCov/circulocov.map-ont.bam"
        # illumina_bam = "CirculoCov/circulocov.sr.bam"

        df_depth = pd.DataFrame(columns = ["contig", "position", "end", "match"])
        df_cov   = pd.DataFrame(columns = ["#rname", "startpos", "endpos"])

        if os.path.exists(nanopore_bam):
            analysis_df_depth, analysis_df_cov = counts(nanopore_bam, genome_dict, "nanopore", args)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   "nanopore")
            if args.all :
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, "nanopore")

        if os.path.exists(illumina_bam):
            analysis_df_depth, analysis_df_cov = counts(illumina_bam, genome_dict, "illumina", args)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   "illumina")
            if args.all:
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, "illumina")

        if os.path.exists(pacbio_bam):
            analysis_df_depth, analysis_df_cov = counts(pacbio_bam, genome_dict, "pacbio", args)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   "pacbio")
            if args.all:
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, "pacbio")

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 4. Extract fastq         #####
        ##### ----- ----- ----- ----- ----- #####

        results_dict = {}
        results_dict['total_length'] = total_length

        if os.path.exists(nanopore_bam):
            results_dict["unmapped_nanopore"] = extract(nanopore_bam, genome_dict,
                                                        "nanopore",
                                                        args,
                                                        temp_dir)

        if os.path.exists(illumina_bam):
            results_dict["unmapped_illumina"] = extract(illumina_bam, genome_dict,
                                                        "illumina",
                                                        args,
                                                        temp_dir)

        if os.path.exists(pacbio_bam):
            results_dict["unmapped_pacbio"]   = extract(pacbio_bam,
                                                        genome_dict,
                                                        "pacbio",
                                                        args,
                                                        temp_dir)

        # used for testing
        #results_dict["unmapped_nanopore"] = 4
        #results_dict["unmapped_illumina"] = 6675

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 5. Graph Coverage        #####
        ##### ----- ----- ----- ----- ----- #####

        if args.all:
            visualize(genome_dict, df_depth, args)

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 6. Create Summary        #####
        ##### ----- ----- ----- ----- ----- #####

        summary(df_cov, genome_dict, results_dict, args)

if __name__ == "__main__":
    main()

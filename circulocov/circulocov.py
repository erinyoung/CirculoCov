#!/usr/bin/env python3

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
import concurrent.futures
import itertools
import logging
import os
import sys
import tempfile
import pandas as pd

#from utils.circular import circular
from utils.genome_prep import genome_prep
from utils.mapping import mapping
from utils.counts import counts
from utils.merge_dataframe import merge_cov_dataframe, merge_depth_dataframe

def main():
    """ Get coverage for draft genomes """

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 0. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    version = '0.1.20230104'

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample', required = False, type = str, help = 'Sample name') # pylint: disable=C0301
    parser.add_argument('-g', '--genome', required = True, type = str, help = 'Genome (draft or complete)') # pylint: disable=C0301
    parser.add_argument('-i', '--illumina', nargs = '+', required = False, type = str, help = 'Input Illumina fastq(s)') # pylint: disable=C0301
    parser.add_argument('-n', '--nanopore', required = False, type = str, help ='Input nanopore fastq') # pylint: disable=C0301
    parser.add_argument('-p', '--pacbio', required = False, type = str, help ='Input pacbio fastq') # pylint: disable=C0301
    parser.add_argument('-d', '--padding', required = False, type = int, help = 'Amount of padding added to circular sequences', default = 10000) # pylint: disable=C0301
    parser.add_argument('-w', '--window', required = False, type = int, help = 'Window size for coverage', default = 500) # pylint: disable=C0301
    parser.add_argument('-o', '--out', required = False, type = str, help = 'directory for results', default = 'CirculoCov') # pylint: disable=C0301
    parser.add_argument('-log', '--loglevel', required = False, type = str, help = 'logging level', default = 'INFO') # pylint: disable=C0301
    parser.add_argument('-t', '--threads', required = False, type = int, help = 'specifies number of threads to use', default=4) # pylint: disable=C0301
    parser.add_argument('-v', '--version', help='print version and exit', action = 'version', version = version) # pylint: disable=C0301
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s - %(message)s',
        datefmt = '%y-%b-%d %H:%M:%S',
        level=args.loglevel.upper())

    if args.sample:
        sample = args.sample
    else:
        sample = "circulocov"

    logging.info(f"Filenames :\t{sample}")
    logging.info('CirculoCov ver :\t'     + str(version))        # pylint: disable=W1201
    logging.info('Final directory :\t' + str(args.out))       # pylint: disable=W1201
    logging.info('Genome file :\t' + str(args.genome))
    logging.info('Window size :\t' + str(args.window))
    logging.info('Num threads :\t' + str(args.threads))
    logging.info('Padding length :\t' + str(args.padding))
    if args.nanopore:
        logging.info('Nanopore file :\t' + str(args.nanopore))
    if args.illumina:
        logging.info('Illumina file(s) :\t' + ', '.join(args.illumina))
    if args.pacbio:
        logging.info('PacBio file :\t' + str(args.pacbio))

    genome   = args.genome
    nanopore = args.nanopore
    illumina = args.illumina
    pacbio   = args.pacbio

    if not nanopore and not illumina and not pacbio:
        logging.fatal('Cannot run without fastq files!')
        sys.exit(1)

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    temp_dir  = tempfile.TemporaryDirectory(dir = args.out) # pylint: disable=R1732
    tmp = temp_dir.name + '/' + sample + '.'
    out = args.out      + '/' + sample + '.'

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 1. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    # Determine if genome is circular
    genome_dict, fasta = genome_prep(genome, tmp, args.padding)
    
    logging.debug("The dictionary for the genome is")
    logging.debug(genome_dict)

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 2. Minimap2              #####
    ##### ----- ----- ----- ----- ----- #####

    # TODO : uncomment these
    # if nanopore:
    #     nanopore_bam = mapping(nanopore, fasta, "map-ont", args.threads, out, tmp)
    # else:
    # nanopore_bam = ""

    # if illumina:
    #     illumina_bam = mapping(illumina, fasta, "sr",      args.threads, out, tmp)
    #    else:
    # illumina_bam = ""

    if pacbio:
        pacbio_bam   = mapping(pacbio,   fasta, "map-pb",  args.threads, out, tmp)
    else:
        pacbio_bam  = ""

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 3. Samtools Counts       #####
    ##### ----- ----- ----- ----- ----- #####

    # TODO : remove these lines (only exists for testing)

    nanopore_bam = "CirculoCov/circulocov.map-ont.bam"
    illumina_bam = "CirculoCov/circulocov.sr.bam"

    df_depth    = pd.DataFrame(columns = ["contig", "position", "end", "match"])
    df_cov      = pd.DataFrame(columns = ["#rname", "startpos", "endpos"])

    nanopore_df_depth = pd.DataFrame()
    illumina_df_depth = pd.DataFrame()
    pacbio_df_depth   = pd.DataFrame()
    nanopore_df_cov   = pd.DataFrame()
    illumina_df_cov   = pd.DataFrame()
    pacbio_df_cov   = pd.DataFrame()

    if os.path.exists(nanopore_bam):
        nanopore_df_depth, nanopore_df_cov = counts(nanopore_bam, genome_dict, args)
        df_depth    = merge_depth_dataframe(df_depth, nanopore_df_depth, "nanopore")
        df_cov      = merge_cov_dataframe(df_cov, nanopore_df_cov, "nanopore")

    exit(0)

    if os.path.exists(illumina_bam):
        illumina_df_depth, illumina_df_cov = counts(illumina_bam, genome_dict, args)
        df_depth       = merge_depth_dataframe(df_cov, illumina_df_depth, "illumina")

    if os.path.exists(pacbio_bam):
        pacbio_df_depth, pacbio_df_cov   = counts(pacbio_bam, genome_dict, args)
        df_depth          = merge_depth_dataframe(df_cov, pacbio_df_depth, "pacbio")


    print("****************")
    print(nanopore_df_depth)
    print(illumina_df_depth)
    print(df_depth)


    ##### ----- ----- ----- ----- ----- #####
    ##### Part 4. Extract fastq         #####
    ##### ----- ----- ----- ----- ----- #####

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 5. Graph Coverage        #####
    ##### ----- ----- ----- ----- ----- #####


    ##### ----- ----- ----- ----- ----- #####
    ##### Part 6. Create Summary        #####
    ##### ----- ----- ----- ----- ----- #####




if __name__ == "__main__":
    main()

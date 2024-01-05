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

# I tried to keep dependencies down...
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

def main():
    """ Get coverage for draft genomes """

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 0. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    version = '0.1.20230104'

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--genome', required = True, type = str, help = 'Genome (draft or complete)') # pylint: disable=C0301
    parser.add_argument('-i', '--illumina', nargs = '+', required = False, type = str, help = 'Input Illumina fastq(s)') # pylint: disable=C0301
    parser.add_argument('-n', '--nanopore', required = False, type = str, help ='Input nanopore fastq') # pylint: disable=C0301
    parser.add_argument('-o', '--out', required = False, type = str, help = 'directory for results', default = 'CirculoCov') # pylint: disable=C0301
    parser.add_argument('-log', '--loglevel', required = False, type = str, help = 'logging level', default = 'INFO') # pylint: disable=C0301
    #parser.add_argument('-t', '--threads', required = False, type = int, help = 'specifies number of threads to use', default=4) # pylint: disable=C0301
    parser.add_argument('-v', '--version', help='print version and exit', action = 'version', version = version) # pylint: disable=C0301
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s - %(message)s',
        datefmt = '%y-%b-%d %H:%M:%S',
        level=args.loglevel.upper())

    logging.info('CirculoCov ver :\t'     + str(version))        # pylint: disable=W1201
    logging.info('Final directory :\t' + str(args.out))       # pylint: disable=W1201
    logging.info('Genome file :\t' + str(args.genome))
    if args.nanopore:
        logging.info('Nanopore file :\t' + str(args.nanopore))
    if args.illumina:
        logging.info('Illumina file(s) :\t' + ', '.join(args.illumina))

    genome   = args.genome
    out      = args.out
    nanopore = args.nanopore
    illumina = args.illumina

    if not nanopore and not illumina:
        logging.fatal('Cannot run without fastq files!')
        sys.exit(1)

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    temp_dir  = tempfile.TemporaryDirectory(dir = args.out) # pylint: disable=R1732
    tmp = temp_dir.name + '/'

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 1. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    # Determine if genome is circular
    genome_dict, fasta = genome_prep(genome, tmp)
    
    logging.debug("The dictionary for the genome is")
    logging.debug(genome_dict)
    
    map_illumina(illumina, fasta)
    #map_nanopore(nanopore, fasta)

    #circular(genome_dict, genome, tmp)

    


if __name__ == "__main__":
    main()

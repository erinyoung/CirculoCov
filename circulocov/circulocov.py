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
import subprocess
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
    ''' Get coverage for draft genomes '''

    ##### ----- ----- ----- ----- ----- #####
    ##### Part 0. Setup                 #####
    ##### ----- ----- ----- ----- ----- #####

    version = '0.1.20240216'

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample',
                        required = False,
                        type = str,
                        help = 'Sample name',
                        default= 'circulocov')
    parser.add_argument('-g', '--genome',
                        required = True,
                        type = str,
                        help = 'Genome (draft or complete)')
    parser.add_argument('-i', '--illumina',
                        nargs = '+',
                        required = False,
                        type = str,
                        help = 'Input illumina fastq(s)')
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
                        help = 'Number of windows for coverage',
                        default = 100)
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

    try:
        result = subprocess.run(['minimap2', '--version'], capture_output=True, text=True, check=True)
        minimap2_ver = result.stdout.strip()
    except subprocess.CalledProcessError as e:
        logging.fatal('Minimap2 not found!')
        logging.fatal(e)
        sys.exit(1)

    # printing everything to the screen
    logging.info(f'Filename prefix :\t{args.sample}')
    logging.info(f'CirculoCov ver :\t{str(version)}')
    logging.info(f'minimap2 ver :\t{str(minimap2_ver)}')
    logging.info(f'Final directory :\t{str(args.out)}')
    logging.info(f'Num threads :\t{str(args.threads)}')
    logging.info(f'Padding length :\t{str(args.padding)}')
    if args.all:
        logging.info('All is set :\tWill create windows and graph coverages')
        logging.info(f'Window number :\t{str(args.window)}')
    logging.info(f'Genome file :\t{str(args.genome)}')
    if args.nanopore:
        logging.info(f'Nanopore file :\t{str(args.nanopore)}')
    if args.illumina:
        logging.info(f"Illumina file(s) :\t{', '.join(args.illumina)}")
    if args.pacbio:
        logging.info(f'PacBio file :\t{str(args.pacbio)}')

    if not args.nanopore and not args.illumina and not args.pacbio:
        logging.fatal('Cannot run without fastq files!')
        sys.exit(1)

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if not os.path.exists(args.out + '/fastq') and args.all :
        os.mkdir(args.out + '/fastq')

    with tempfile.TemporaryDirectory(dir = args.out) as temp_dir:

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 1. Setup                 #####
        ##### ----- ----- ----- ----- ----- #####
        
        logging.info('Setting up genome file')
        genome_dict, fasta, total_length = genome_prep(args, temp_dir)

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 2. Minimap2              #####
        ##### ----- ----- ----- ----- ----- #####

        logging.info('Mapping reads to reference')

        if args.nanopore:
            logging.info('Mapping nanopore reads to reference')
            nanopore_bam = mapping(args.nanopore, fasta, 'map-ont', args, temp_dir)
        else:
            nanopore_bam = ''

        if args.illumina:
            logging.info('Mapping illumina reads to reference')
            illumina_bam = mapping(args.illumina, fasta, 'sr', args, temp_dir)
        else:
            illumina_bam = ''

        if args.pacbio:
            logging.info('Mapping pacbio reads to reference')
            pacbio_bam   = mapping(args.pacbio, fasta, 'map-pb', args, temp_dir)
        else:
            pacbio_bam   = ''

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 3. Samtools Counts       #####
        ##### ----- ----- ----- ----- ----- #####

        logging.info('Getting coverage and depth')

        df_depth = pd.DataFrame(columns = ['contig', 'pos'])
        df_cov   = pd.DataFrame(columns = ['#rname', 'startpos', 'endpos'])

        if os.path.exists(nanopore_bam):
            logging.info('Getting coverage and depth for nanopore reads')
            analysis_df_depth, analysis_df_cov = counts(nanopore_bam, genome_dict, 'nanopore', args, temp_dir)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   'nanopore')
            if args.all :
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, 'nanopore')

        if os.path.exists(illumina_bam):
            logging.info('Getting coverage and depth for illumina reads')
            analysis_df_depth, analysis_df_cov = counts(illumina_bam, genome_dict, 'illumina', args, temp_dir)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   'illumina')
            if args.all:
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, 'illumina')

        if os.path.exists(pacbio_bam):
            logging.info('Getting coverage and depth for PacBio reads')
            analysis_df_depth, analysis_df_cov = counts(pacbio_bam, genome_dict, 'pacbio', args, temp_dir)
            df_cov   = merge_cov_dataframe(  df_cov,   analysis_df_cov,   'pacbio')
            if args.all:
                df_depth = merge_depth_dataframe(df_depth, analysis_df_depth, 'pacbio')
        
        df_depth = df_depth.infer_objects(copy=False).fillna(0)
        df_depth = df_depth.sort_values(by=['contig', 'pos']).reset_index(drop=True)
        df_depth.to_csv(args.out + '/depth.txt', index=False, sep = '\t')
        
        df_cov = df_cov.sort_values(by=['endpos'], ascending=[False]).reset_index(drop=True)
        df_cov.to_csv(  args.out + '/cov.txt',   index=False, sep = '\t')

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 4. Extract fastq         #####
        ##### ----- ----- ----- ----- ----- #####
                
        logging.info('Extracting fastq files')

        results_dict = {}
        results_dict['total_length'] = total_length

        if os.path.exists(nanopore_bam):
            logging.info('Extracting nanopore fastq files')
            results_dict['unmapped_nanopore'] = extract(nanopore_bam, genome_dict,
                                                        'nanopore',
                                                        args,
                                                        temp_dir)

        if os.path.exists(illumina_bam):
            logging.info('Extracting illumina fastq files')
            results_dict['unmapped_illumina'] = extract(illumina_bam, genome_dict,
                                                        'illumina',
                                                        args,
                                                        temp_dir)

        if os.path.exists(pacbio_bam):
            logging.info('Extracting pacbio fastq files')
            results_dict['unmapped_pacbio']   = extract(pacbio_bam,
                                                        genome_dict,
                                                        'pacbio',
                                                        args,
                                                        temp_dir)

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 5. Graph Coverage        #####
        ##### ----- ----- ----- ----- ----- #####

        df_depth = pd.read_table(args.out + '/depth.txt')

        if args.all:
            logging.info('Graphing coverage')

            visualize(genome_dict, df_depth, args)

        ##### ----- ----- ----- ----- ----- #####
        ##### Part 6. Create Summary        #####
        ##### ----- ----- ----- ----- ----- #####
            
        logging.info('Creating summary')
        summary(df_cov, genome_dict, results_dict, args)

if __name__ == '__main__':
    main()

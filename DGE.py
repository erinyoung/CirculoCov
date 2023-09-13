#!/usr/bin/env python3

'''
Author: Erin Young

Description:

This script takes nanopore and/or illumina reads and maps them to a fasta file with minimap2.

This generates a bam file which is then used determine the coverage of each fasta file.

Two files are generated for this coverage:
- <outdir>/coverage.txt is a tab-delimited file with each contig and the values of coverage over windows.
- <outdir>/coverage.png is a figure of the largest 10 contigs and their nanopore/illumina coverage

Tis masqurading as a real python program.

EXAMPLE:
dgc -r input.fasta -i illumina_R1.fastq illumina_R2.fastq -n nanopore.fastq -o outdir
'''

# I tried to keep dependencies down...
import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import mappy as mp
import subprocess
from pybedtools import BedTool
from shutil import which


#logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%y-%b-%d %H:%M:%S', level=logging.INFO)
version = '0.0.20230912'

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%y-%b-%d %H:%M:%S', level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--illumina',  type=str, nargs = '+', default = '', help = 'illumina reads in fastq or fastq.gz format')
parser.add_argument('-n', '--nanopore',  type = str, default = '', help ='nanopore reads in fastq or fastq.gz format')
parser.add_argument('-r', '--reference', type = str, required = True, help = 'fasta file to map/align reads to')
parser.add_argument('-o', '--out',       type=str, help = 'directory for results', default='dge')
parser.add_argument('-w', '--window',    type = int, help = 'window size for coverage', default = 500)
parser.add_argument('-t', '--threads',   type=int, help = 'specifies number of threads to use', default=4)
parser.add_argument('-v', '--version',   help = 'print version and exit', action = 'version', version = '%(prog)s ' + version)
args   = parser.parse_args()

def align( reference, reads, threads, out, type):
    filename = out + '/tmp.' + type + '.sam'
    with open(filename, "w") as sam:
        run = subprocess.run(['minimap2', '-ax', type, '-t', str(threads), reference] + reads, stdout = sam)
        logging.debug('The commandline is {}'.format(run.args))
    pysam.sort('-o', out + '/tmp.' + type + '.sorted.sam', '-@', str(threads), filename)

def check_tool(tool):
    path = which(tool)
    logging.debug('The path for each tool is ' + str(path))
    if path == None :
        logging.info('Missing : ' + tool)
        found = 1
    else:
        logging.info('Found : ' + tool + ' at ' + str(path))
        found = 0
    return(found)

def get_reference_length(fasta, window, out):
    pysam.faidx(fasta)
    df_fai = pd.read_table(fasta + '.fai', usecols=[0,1], header = None)
    df_fai.columns = ['NAME', 'LENGTH']
    logging.debug(df_fai)

    df_ref = pd.DataFrame(columns = ['#name', 'start', 'end'])
    for ind in df_fai.index:
        contig        = df_fai.at[df_fai.index[ind], 'NAME']
        contig_length = df_fai.at[df_fai.index[ind], 'LENGTH']
        logging.debug('Contig ' + str(contig) + ' has the length of ' + str(contig_length))
        start = 0
        while start < contig_length:
            end = start + window - 1
            df_ref.loc[len(df_ref)] = {'#name': contig, 'start': start, 'end': end}
            start = start + window

    bed = BedTool.from_dataframe(df_ref)
    df_ref.to_csv(out + '/tmp.bed', index=False, sep='\t')
    return bed

def get_window_coverage(bam, out, bed):
    #bed = out + '/tmp.bed'
    #pysam.depth('-b', bed, '-o', bam + '.cov', bam)
    #bed.BedTool.multi_bam_coverage(bam)
    bedtools_bed = BedTool('dge/tmp.bed')
    print(bedtools_bed)
    bedtools_bed.multi_bam_coverage(bams='dge/tmp.map-ont.sorted.sam')
    #bedtools_bam = Bedtool('dge/tmp.map-ont.sorted.sam')
    
    #multi_bam_coverage(bams=bedtools_bam, bed=bedtools_bed)
    #print(cov)

if __name__ == "__main__":
    # check that input files were specified (does not check if they exist)
    if not args.illumina and not args.nanopore:
        logging.error('No illumina or nanopore reads specified. Exiting')
        exit(1)

    # make the final directory
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    logging.info('DGC version :\t\t'          + str(version))
    logging.info('Number of threads :\t'      + str(args.threads))
    logging.info('Final directory :\t\t'      + str(args.out))
    if args.illumina: logging.info('Input illumina file(s) :\t' + ', '.join(args.illumina))
    if args.nanopore: logging.info('Input nanopore file :\t'    + str(args.nanopore))
    logging.info('Input fasta file : \t'      + str(args.reference))

    reference = args.reference
    illumina  = list(args.illumina)
    nanopore  = [ args.nanopore]
    window    = args.window
    threads   = args.threads
    out       = args.out
    
    logging.debug(illumina)
    logging.debug(nanopore)
    
    illuminasam = ''
    nanoporesam = ''

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### checking shell tools                #####
    ##### ----- ----- ----- ----- ----- ----- #####
    
    tools = ['minimap2']
    find_tool = 0
    for tool in tools:
        check = check_tool(tool)
        if check == 1:
            find_tool = 1
        logging.debug('The check boolean is ' + str(find_tool))

    if find_tool > 0 :
        logging.error('Dependencies were not found!')
        exit(1)

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### aligning reads to fasta             #####
    ##### ----- ----- ----- ----- ----- ----- #####

#    if illumina:
#        logging.info('Starting alignment for Illumina reads')
#        align( reference, illumina, threads, out, 'sr')
#        illuminasam = out + '/tmp.sr.sam'

#    if nanopore:
#        logging.info('Starting alignment for ONT reads')
#        align( reference, nanopore, threads, out, 'map-ont')
#        nanoporesam = out + '/tmp.map-ont.sam'

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### getting coverage                    #####
    ##### ----- ----- ----- ----- ----- ----- #####

    bed = get_reference_length(reference, window, out)
    #out + '/tmp.bed'

    # step 3 : get coverage for each window
    get_window_coverage(out + '/tmp.map-ont.sorted.sam', out, bed)


    # step 4 : get average coverage for each contig
    # step 5 : get average coverage for all contigs

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### creating figures                    #####
    ##### ----- ----- ----- ----- ----- ----- #####

    # from pandas dataframe graph each contig with something
    # graph coverage depth with track for illumina and track for nanopore
    # set line for average for contig
    # set line for average for all

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### removing tmp files                  #####
    ##### ----- ----- ----- ----- ----- ----- #####

    # os.remove(bam0)


    logging.info('DGC is complete! Final files are in args.out')
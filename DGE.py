#!/usr/bin/env python3

'''
Author: Erin Young

Description:

This script takes nanopore and/or illumina reads and maps them to a fasta file or gfa file with minimap2.

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
import concurrent.futures
import itertools
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import mappy as mp
import subprocess
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
    sorted_filename = out + '/tmp.' + type + '.sorted.bam'
    with open(filename, "w") as sam:
        run = subprocess.run(['minimap2', '-ax', type, '-t', str(threads), reference] + reads, stdout = sam)
        logging.debug('The commandline is {}'.format(run.args))
    pysam.sort('-o', sorted_filename, '-@', str(threads), filename)
    pysam.index(sorted_filename)
    return(sorted_filename)

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

def convert_to_fasta(gfa):
    print(gfa)

def get_coverage(df, bam, region):
    row = int(region.split('-')[0].split(':')[1]) 
    logging.debug('Getting coverage for ' + bam + ' and ' + region + ' on row ' + str(row))
    cov=float(pysam.coverage('--no-header', bam, '-r', region).split()[6])
    logging.debug('The coverage is ' + str(cov))
    rowindex = df.index[df['region'] == row ]
    df.loc[rowindex, [bam]] = cov

def get_bam_coverage(df_bam, bam):
    logging.debug('Getting coverage for entire ' + bam)
    cov=float(pysam.coverage('--no-header', bam).split()[6])
    
    total    = float(pysam.view('-c', bam).strip())    
    unmapped = float(pysam.view('-c', '-f', '4', bam).strip())
    mapped   = float(pysam.view('-c', '-F', '4', bam).strip())

    rowindex = df_bam.index[df_bam['bam'] == bam ]
    df_bam.loc[rowindex, ['cov']]      = cov
    df_bam.loc[rowindex, ['unmapped']] = unmapped
    df_bam.loc[rowindex, ['mapped']]   = mapped
    df_bam.loc[rowindex, ['total']]    = total

def get_contig_coverage(df_contig, bam, region):
    logging.debug('Getting coverage for ' + bam + ' and ' + region)
    cov=float(pysam.coverage('--no-header', bam, '-r', region).split()[6])
    logging.debug('The coverage is ' + str(cov))
    rowindex = df_contig.index[df_contig['region'] == region ]
    df_contig.loc[rowindex, [bam]] = cov

def get_reference_length(fasta, window, out):
    pysam.faidx(fasta)
    df_fai = pd.read_table(fasta + '.fai', usecols=[0,1], header = None)
    df_fai.columns = ['NAME', 'LENGTH']
    logging.debug(df_fai)

    ranges  = []
    contigs = []
    starts  = []

    for ind in df_fai.index:
        contig        = df_fai.at[df_fai.index[ind], 'NAME']
        contigs.append(contig)
        contig_length = df_fai.at[df_fai.index[ind], 'LENGTH']
        logging.debug('Contig ' + str(contig) + ' has the length of ' + str(contig_length))
        start = 0
        while start < contig_length:
            end = start + window - 1
            new_range = contig + ':' + str(start) + '-' + str(end)
            ranges.append( new_range )
            starts.append( start )
            start = start + window

    return ranges, contigs, starts

if __name__ == "__main__":
    # check that input files were specified (does not check if they exist)
    if not args.illumina and not args.nanopore:
        logging.error('No illumina or nanopore reads specified. Exiting')
        exit(1)

    # make the final directory
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    logging.info('DGC version :\t\t'     + str(version))
    logging.info('Final directory :\t\t' + str(args.out))
    logging.info('Number of threads :\t' + str(args.threads))
    logging.info('Window size : \t\t'    + str(args.window))
    logging.info('Input fasta file : \t' + str(args.reference))
    if args.illumina: logging.info('Input illumina file(s) :\t' + ', '.join(args.illumina))
    if args.nanopore: logging.info('Input nanopore file :\t'    + str(args.nanopore))

    reference = args.reference
    illumina  = list(args.illumina)
    nanopore  = [args.nanopore]
    window    = args.window
    threads   = args.threads
    out       = args.out
    
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
    ##### converting gfa to fasta             #####
    ##### ----- ----- ----- ----- ----- ----- #####

    if 'gfa' in reference:
        logging.info('Conferting gfa file ' + reference + ' to fasta file')
        reference, extra_reference = convert_to_fasta(gfa)
        regions, contigs, starts = get_reference_length(extra_reference, window, out)
    else:
        regions, contigs, starts = get_reference_length(reference, window, out)

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### aligning reads to fasta             #####
    ##### ----- ----- ----- ----- ----- ----- #####

    bams = []
    columns = ['region']

    if illumina:
        logging.info('Starting alignment for Illumina reads')
        bam = align( reference, illumina, threads, out, 'sr')
        bams.append(bam)
        columns.append('illumina')

    if nanopore:
        logging.info('Starting alignment for ONT reads')
        bam = align( reference, nanopore, threads, out, 'map-ont')
        bams.append(bam)
        columns.append('nanopore')

    logging.debug('The bams :')
    logging.debug(bams)

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### getting coverage                    #####
    ##### ----- ----- ----- ----- ----- ----- #####

    # bams = ['dge/tmp.sr.sorted.bam', 'dge/tmp.map-ont.sorted.bam'] # print

    logging.info('Getting average coverage for each window')
    df_windows           = pd.DataFrame([])
    df_windows['region'] = starts
    pool = concurrent.futures.ThreadPoolExecutor( max_workers = threads )
    for input in list(itertools.product( bams, regions )):
        pool.submit(get_coverage, df_windows, input[0], input[1])
        # get_coverage(df_windows, input[0], input[1])

    logging.info('Getting average coverage for each contig')
    df_contigs = pd.DataFrame([])
    df_contigs['region'] = contigs
    for input in list(itertools.product( bams, contigs )):
        # get_contig_coverage(df_contigs, input[0], input[1])
        pool.submit(get_contig_coverage, df_contigs, input[0], input[1]) 

    logging.info('Getting average coverage')
    df_bams = pd.DataFrame([])
    df_bams['bam'] = bams
    for bam in bams:
        # get_bam_coverage(df_bams, bam)
        pool.submit(get_bam_coverage, df_bams, bam)

    pool.shutdown(wait=True)

    columns = ['region', 'illumina', 'nanopore']
    df_windows.columns = columns
    df_contigs.columns = columns
    df_bams['type'] = df_bams['bam'].str.replace(out + '/', '', regex = False).replace('tmp.sr.sorted.bam', 'illumina', regex = False).replace('tmp.map-ont.sorted.bam', 'nanopore', regex = False)
    logging.debug(df_windows)
    logging.debug(df_contigs)
    logging.debug(df_bams)

    # writing results to a file
    df_windows.to_csv(out + '/windows_depth.csv', index=False)
    df_contigs.to_csv(out + '/contigs_depth.csv', index=False)
    df_bams.to_csv(   out + '/average_depth.csv', index=False)

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### removing tmp files                  #####
    ##### ----- ----- ----- ----- ----- ----- #####

    tmp_files = [ 'tmp.map-ont.sam', 'tmp.map-ont.sorted.bam', 'tmp.map-ont.sorted.bam.bai', 'tmp.sr.sam', 'tmp.sr.sorted.bam', 'tmp.sr.sorted.bam.bai' ]
    for file in tmp_files:
        os.remove(out + '/' + file)

    ##### ----- ----- ----- ----- ----- ----- #####
    ##### creating figures                    #####
    ##### ----- ----- ----- ----- ----- ----- #####

    # from pandas dataframe graph each contig with something
    # graph coverage depth with track for illumina and track for nanopore
    # set line for average for contig
    # set line for average for all

    # getting rid of column with all the bam names
    #df.drop('bam', axis=1, inplace=True)

    #boxplot = df.boxplot(fontsize=5, rot=90, figsize=(15,8), grid=False)
    #boxplot.plot()
    #plt.title('Primer Assessment')
    #boxplot.set_ylabel('meandepth')
    #boxplot.set_xlabel('amplicon name')
    #boxplot.figure.savefig(args.out + '/amplicon_depth.png', bbox_inches='tight')
    #plt.close()

    logging.info('DGC is complete! Final files are in ' + out)
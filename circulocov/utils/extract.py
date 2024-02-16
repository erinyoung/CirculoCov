#!/usr/bin/env python
# pylint: disable=logging-fstring-interpolation

''' Extract reads for each contig '''

import logging
import pysam

def extract(bam, genome_dict, analysis, args, temp_dir):
    ''' Extract reads for each contig '''

    if args.all :
        # Checking if single end or paired end
        single_check = int(pysam.view('-c', '-f',  '1', bam))

        # NOTE : pysam sort has it's own multithreading, so don't use concurrent
        for contig in genome_dict.keys():
            logging.debug(f'Extracting reads for {contig} from {analysis} bam file')

            # bam file names
            out_bam    = temp_dir + '/' + args.sample + '_' + contig + '_' + analysis + '.bam'
            sorted_bam = temp_dir + '/' + args.sample + '_' + contig + '_' + analysis + '.bam'

            # fastq file names
            fastq_base  = args.out + '/fastq/' + args.sample + '_' + contig + '_' + analysis

            # getting mapped reads for region
            pysam.view('-b', '-F', '4', bam, '-o', out_bam, contig, catch_stdout=False)
            pysam.sort('-o', sorted_bam, '-@', str(args.threads), out_bam, catch_stdout=False)
            pysam.index(sorted_bam)
            if single_check == 0 :
                pysam.fastq('-@',
                            str(args.threads),
                            sorted_bam,
                            '-0',
                            fastq_base + '.fastq.gz',
                            catch_stdout=False)
            else:
                pysam.fastq('-@',
                            str(args.threads),
                            sorted_bam,
                            '-1', fastq_base + '_R1.fastq.gz',
                            '-2', fastq_base + '_R2.fastq.gz',
                            '-s', fastq_base + '_singletons.fastq.gz', 
                            catch_stdout=False)

        # get unmapped reads as well
        logging.debug(f'Extracting unmapped reads from {analysis} bam file')
        unmapped_bam        = temp_dir + '/' + args.sample + '_unmapped_' + analysis + '.bam'
        sorted_unmapped_bam = temp_dir + '/' + args.sample + '_unmapped_' + analysis + '.bam'
        unmapped_base       = args.out + '/fastq/' + args.sample + '_unmapped_' + analysis

        pysam.view('-b', '-f', '4', bam, '-o', unmapped_bam, catch_stdout=False)
        pysam.sort('-o', sorted_unmapped_bam, '-@', str(args.threads), unmapped_bam, catch_stdout=False)
        pysam.index(sorted_unmapped_bam)
        if single_check == 0 :
            pysam.fastq('-@',
                        str(args.threads),
                        sorted_unmapped_bam,
                        '-0', unmapped_base + '.fastq.gz',
                        catch_stdout=False)
        else:
            pysam.fastq('-@',
                        str(args.threads),
                        sorted_unmapped_bam,
                        '-1', unmapped_base + '_R1.fastq.gz',
                        '-2', unmapped_base + '_R2.fastq.gz',
                        '-s', unmapped_base + '_singletons.fastq.gz',
                        catch_stdout=False)

    num_unmapped = int(pysam.view('-c', '-f',  '4', bam))
    logging.info(f'There are {str(num_unmapped)} unmapped {analysis} reads.')

    return num_unmapped

#!/usr/bin/env python

""" Add sequences to the end for circular contigs """

import logging
import os
import pysam
import shutil
from Bio import SeqIO

def circular(genome_dict, genome, out):
    """ Add sequences to end if circular """

    #final_genome = out + "/" + os.path.basename(genome)
    #shutil.copyfile(genome, final_genome)

    for contig in genome_dict.keys():
        print(contig)
        print(genome_dict[contig])
        pysam.faidx(genome, contig, "-o", out + contig + '.1.fasta')
        #if "circular=true" in genome_dict[contig]["Description"]:
        #    print(genome_dict[contig])
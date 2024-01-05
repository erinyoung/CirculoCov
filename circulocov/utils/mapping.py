#!/usr/bin/env python

""" Add sequences to the end for circular contigs """

import logging
import os
import pysam
import shutil
from Bio import SeqIO

def minimap(genome_dict, genome, out):
    """ Minimap2 FTW """

    subprocess.run(minimap2, '-d', , target_file  + threads )
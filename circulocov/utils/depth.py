#!/usr/bin/env python

''' Gets depth for bam '''

import pysam

def depth(bam, contig, tmp):
    ''' Gets coverage for bam for contig '''

    depth_file = tmp + '/' + contig + '.txt'

    pysam.depth(bam, '-r', contig, '-o', depth_file )

    return depth_file

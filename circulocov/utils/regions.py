#!/usr/bin/env python

""" Creating region list """

import logging

def regions(genome_dict, args):
    """ Creating region list """
    regs = []

    for contig in genome_dict.keys():
        for i in range(1, genome_dict[contig]['length'], args.window):
            end = i + args.window - 1
            if end > genome_dict[contig]['length']:
                end = genome_dict[contig]['length']
            regs.append(contig + ":" + str(i) + "-" + str(end))

    logging.debug("The regions are")
    logging.debug(regs)

    return regs

def padded_regions(genome_dict, args):
    """ Creating padded region list """

    regs = []
    for contig in genome_dict.keys():
        if genome_dict[contig]['length'] < args.padding :
            padded_length = genome_dict[contig]['length'] + genome_dict[contig]['length']
        else:
            padded_length = genome_dict[contig]['length'] + args.padding

        for i in range(genome_dict[contig]['length'] + 1, padded_length, args.window):
            end = min(i + args.window - 1, padded_length)

            regs.append(contig + ":" + str(i) + "-" + str(end))

    logging.debug("The padded regions are:")
    logging.debug(regs)

    return regs

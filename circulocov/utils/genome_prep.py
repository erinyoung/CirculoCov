#!/usr/bin/env python
# pylint: disable=logging-fstring-interpolation

''' Prepare genome for assembly '''

import logging
import os
from Bio import SeqIO
#from Bio.SeqUtils import GC

def genome_prep(args, out):
    ''' Turn fasta file into dict '''

    prepped_fasta = out + '/' + os.path.basename(args.genome)

    # to be returned to script
    genome_dict = {}

    # for appending together
    records = []
    circ_count = 0
    total_length = 0
    for record in SeqIO.parse(args.genome,'fasta'):
        genome_dict[record.id] = {}
        genome_dict[record.id]['length'] = len(record.seq)
        genome_dict[record.id]['name'] = str(record.id)
        #genome_dict[record.id]['GC'] = GC(record.seq)
        total_length = total_length + len(record.seq)

        # checking if circular
        circular_trues = ['circular=true', 'circ=true', 'circular=t', 'circ=t', 'complete sequence']
        if any(x in record.description.lower() for x in circular_trues):
            genome_dict[record.id]['circ'] = True
            circ_count += 1
            if len(record.seq) < args.padding:
                record.seq = record.seq + record.seq
            else:
                record.seq = record.seq + record.seq[:args.padding + 1]
        else:
            genome_dict[record.id]['circ'] = False

        records.append(record)

    # creating file for mapping
    SeqIO.write(records, prepped_fasta, 'fasta')

    logging.info(f'There were {str(len(genome_dict.keys()))} sequences found in {args.genome}')
    logging.info(f'There were {str(circ_count)} circular sequences')

    logging.debug('The dictionary for the genome is')
    logging.debug(genome_dict)

    return genome_dict, prepped_fasta, total_length

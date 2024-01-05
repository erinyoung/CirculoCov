#!/usr/bin/env python

""" Prepare genome for assembly """

import logging
import os
from Bio import SeqIO

def genome_prep(fasta, out):
    """ Turn fasta file into dict """

    prepped_fasta = out + os.path.basename(fasta)

    # to be returned to script
    genome_dict = {}

    # for appending together
    records = []
    circ_count = 0
    for record in SeqIO.parse(fasta,"fasta"):
        genome_dict[record.id] = {}
        genome_dict[record.id]["length"] = len(record.seq)
        
        # checking if circular
        if any(x in record.description.lower() for x in ["circular=true", "circ=true", "circular=t", "circ=t"]):
            genome_dict[record.id]["circ"] = True
            circ_count += 1
            if len(record.seq) < 10000:
                record.seq = record.seq + record.seq
            else:
                record.seq = record.seq + record.seq[:10001]
        else:
            genome_dict[record.id]["circ"] = False

        records.append(record)

    # creating file for mapping
    SeqIO.write(records, prepped_fasta, "fasta")

    logging.info("There were " + str(len(genome_dict.keys())) + " sequences found in " + fasta)
    logging.info("There were " + str(circ_count) + " circular sequences")

    return genome_dict, prepped_fasta
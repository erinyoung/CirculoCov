#!/usr/bin/env python

''' Get 10 longest contigs '''

def top_ten (genome_dict):

    ''' Get 10 longest contigs '''

    lengths = []
    for key in genome_dict.keys():
        lengths.append(genome_dict[key]['length'])

    if len(lengths) > 9:
        length_threshold = sorted(lengths, reverse=True)[9]
    else:
        length_threshold = 0

    return length_threshold

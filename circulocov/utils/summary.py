#!/usr/bin/env python

''' Creating summary file '''

from .average_depth import average_depth
from .overall_row import overall_row

def summary(df_cov, genome_dict, results_dict, args):
    ''' Creating summary file '''
    # count if single end or paired end

    results_dict['meandepth'] = average_depth(df_cov)

    df = overall_row(df_cov, results_dict, args)

    df['circ'] = 'X'
    for contig in genome_dict.keys():
        df.loc[df['#rname'] == contig, 'circ'] = genome_dict[contig]['circ']
        df.loc[df['#rname'] == contig, 'endpos'] = genome_dict[contig]['length']

    # now to get the dataframe pretty
    first_column = df.pop('circ')
    df.insert(0, 'circ', first_column)

    first_column = df.pop('sample')
    df.insert(0, 'sample', first_column)

    df = df.sort_values(by=['endpos', '#rname'], ascending= [False, True], ignore_index=True)
    df = df.rename(columns={'#rname': 'contigs', 'endpos': 'length'})

    df.to_csv(args.out + '/overall_summary.txt', index=False, sep = '\t')

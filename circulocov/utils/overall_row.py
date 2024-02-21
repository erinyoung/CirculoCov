#!/usr/bin/env python

''' Creating overall summary '''

import pandas as pd

def overall_row(df_cov, results_dict, args):
    ''' Creating an overall summary '''
    # count if single end or paired end

    df = df_cov.copy()

    # dropping columns
    if 'match' in df.columns:
        df = df.drop('match', axis=1)

    if 'weighted' in df.columns:
        df = df.drop('weighted', axis=1)

    if 'startpos' in df.columns:
        df = df.drop('startpos', axis=1)

    for analysis in ['nanopore', 'illumina', 'pacbio']:
        for header in ['meanbaseq', 'meanmapq']:
            if analysis + '_' + header in df.columns:
                df = df.drop(analysis + '_' + header, axis=1)

    orig_len = len(df)

    df.loc[orig_len] = pd.Series()

    df.loc[orig_len,'#rname'] = 'all'
    df.loc[orig_len,'endpos'] = results_dict['total_length']

    after_len = len(df)
    df.loc[after_len] = pd.Series()
    df.loc[after_len,'#rname'] = 'missing'
    df.loc[after_len,'endpos'] = 1

    df['sample'] = args.sample

    for analysis in ['nanopore', 'illumina', 'pacbio']:
        if analysis + '_numreads' in df.columns:
            df.loc[after_len, analysis + '_numreads'] = results_dict['unmapped_' + analysis ]
            df[analysis + '_numreads'] = pd.to_numeric(df[analysis + '_numreads'])
            summation = df[analysis + '_numreads'].sum()
            df.loc[orig_len, analysis + '_numreads'] = summation + results_dict['unmapped_' + analysis] # pylint disable=C0301
            df[analysis + '_numreads'] = df[analysis + '_numreads'].astype(int)

        for header in ['covbases']:
            if analysis + '_' + header in df.columns:
                df.loc[after_len, analysis + '_' + header ] = 0
                df[analysis + '_' + header] = pd.to_numeric(df[analysis + '_' + header])
                summation = df[analysis + '_' + header].sum()
                df.loc[orig_len, analysis + '_' + header] = summation
                df[analysis + '_' + header] = df[analysis + '_' + header].astype(int)

        for header in ['meandepth', 'coverage']:
            if analysis + '_' + header in results_dict['meandepth'].keys():
                df.loc[after_len, analysis + '_' + header ] = 0
                df.loc[orig_len, analysis + '_' + header ] = results_dict['meandepth'][analysis + '_' + header] # pylint disable=C0301
                df[analysis + '_' + header ] = df[analysis + '_' + header ].round(2)

    return df

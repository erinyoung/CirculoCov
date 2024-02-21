#!/usr/bin/env python

''' Getting the average coverage '''

import logging
import pandas as pd

def weighted_meandepth(df, analysis, header):
    ''' Gets weighted average '''
    df[analysis + '_' + header] = pd.to_numeric(df[analysis + '_' + header], downcast='integer')
    df['weighted'] = df['endpos'] * df[analysis + '_' + header]
    weighted       = df['weighted'].sum()
    total_length   = df['endpos'].sum()
    weighted_ave   = round(weighted / total_length, 2)
    df             = df.drop('weighted', axis=1)
    return weighted_ave

def average_depth(df):
    ''' Gets average depth '''

    depth_dict = {}

    if 'nanopore_meandepth' in df.columns:
        depth_dict['nanopore_meandepth'] = weighted_meandepth(df, 'nanopore', 'meandepth')
        depth_dict['nanopore_coverage']  = weighted_meandepth(df, 'nanopore', 'coverage')

    if 'illumina_meandepth' in df.columns:
        depth_dict['illumina_meandepth'] = weighted_meandepth(df, 'illumina', 'meandepth')
        depth_dict['illumina_coverage']  = weighted_meandepth(df, 'illumina', 'coverage')

    if 'pacbio_meandepth' in df.columns:
        depth_dict['pacbio_meandepth'] = weighted_meandepth(df, 'pacbio', 'meandepth')
        depth_dict['pacbio_coverage']  = weighted_meandepth(df, 'pacbio', 'coverage')

    logging.debug('depth dict:')
    logging.debug(depth_dict)

    return depth_dict

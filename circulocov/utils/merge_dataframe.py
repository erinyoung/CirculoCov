#!/usr/bin/env python
# pylint: disable=logging-fstring-interpolation

''' Making mergine pretty '''

import logging
import pandas as pd

def merge_depth_dataframe(df, analysis_df, analysis):
    ''' Adds prefixes and removes unneeded columns after merging '''

    def fill_contig_column(row):
        if pd.isna(row['contig']):
            return row[analysis + '_contig']
        else:
            return row['contig']

    def fill_pos_column(row):
        if pd.isna(row['pos']):
            return row[analysis + '_pos']
        else:
            return row['pos']

    analysis_df['merge'] = analysis_df['contig'].astype(str) + '_' + analysis_df['pos'].astype(str)
    analysis_df = analysis_df.add_prefix(analysis + '_')

    copy_df = df.copy()
    copy_df['merge'] = copy_df['contig'].astype(str) + '_' + copy_df['pos'].astype(str)

    full_df           = pd.merge(copy_df, analysis_df, left_on='merge', right_on= analysis + '_merge', how='outer')
    full_df['contig'] = full_df.apply(fill_contig_column, axis=1)
    full_df['pos']    = full_df.apply(fill_pos_column,    axis=1)
    full_df           = full_df.drop(['merge', analysis + '_merge', analysis + '_contig', analysis + '_pos'], axis=1)
    full_df['pos']    = pd.to_numeric(full_df['pos'], downcast='integer')

    logging.info(f'Mean depths from {analysis} have been merged in')
    return full_df

def merge_cov_dataframe(df, analysis_df, analysis):
    ''' Adds prefixes and removes unneeded columns after merging '''

    analysis_df['match'] = analysis_df['#rname'] + '-' + analysis_df['startpos'].astype(str) + '-' + analysis_df['endpos'].astype(str) # pylint: disable=C0301
    analysis_df = analysis_df.add_prefix(analysis + '_')

    df['match'] = df['#rname'] + '-' + df['startpos'].astype(str) + '-' + df['endpos'].astype(str)
    df = pd.merge(df, analysis_df, left_on='match', right_on = analysis + '_match', how='outer')

    for header in ['#rname', 'startpos', 'endpos', 'match']:
        df.loc[df[header].isnull(),  header ] = df[analysis + '_' + header ]
        df = df.drop(analysis + '_' + header, axis=1)

    logging.info(f'Coverages from {analysis} have been merged in')
    return df

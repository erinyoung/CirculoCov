#!/usr/bin/env python

''' Get dataframe of coverage '''

import pandas as pd

def create_depth_dataframe(df, genome_dict, args):
    ''' Creating dataframe of coverages '''

    df['contig'] = df['contig'].astype(str)
    df['pos']    = pd.to_numeric(df['pos'],   downcast='integer')
    df['depth']  = pd.to_numeric(df['depth'], downcast='integer')

    df_window_depth = pd.DataFrame(columns = ['contig', 'pos', 'depth'])

    for contig in genome_dict.keys():

        divisor = round(genome_dict[contig]['length'] / args.window)

        # getting the 'divisor' points on the graph plus the near beginning and near end       
        df_filtered = df[(df['contig'] == contig) & (df['pos'] <= genome_dict[contig]['length']) & ((df['pos'] % divisor == 0) | ( df['pos'] == 1 ) | ( df['pos'] == genome_dict[contig]['length'] - 1 ) )].copy()
        if 1 not in df_filtered['pos'].values:
            df_filtered.loc[len(df_filtered.index)] = [contig, 1, 0]  
        if genome_dict[contig]['length'] - 1 not in df_filtered['pos'].values:
            df_filtered.loc[len(df_filtered.index)] = [contig, genome_dict[contig]['length'] - 1, 0]

        # adding the padded coverage
        df_filtered_padded = df[(df['contig'] == contig) & (df['pos'] > int(genome_dict[contig]['length']))].copy()
        if not df_filtered_padded.empty :
            df_filtered_padded['new_pos'] = df_filtered_padded['pos'] - int(genome_dict[contig]['length'])
            df_filtered_padded = df_filtered_padded[(df_filtered_padded['new_pos'] % divisor == 0 ) | (df_filtered_padded['new_pos'] == 1 ) | (df_filtered_padded['new_pos'] == genome_dict[contig]['length'] - 1 ) ].copy()
            df_filtered_padded = df_filtered_padded.drop(['pos', 'contig'], axis=1)
            df_filtered_padded = df_filtered_padded.rename(columns={'depth': 'padded_depth', 'new_pos': 'pos'})

            df_filtered = pd.merge(df_filtered, df_filtered_padded, left_on='pos', right_on='pos', how='left')
            df_filtered = df_filtered.fillna(0)
            df_filtered['new_depth'] = df_filtered['padded_depth'] + df_filtered['depth']
            df_filtered = df_filtered.drop(['padded_depth', 'depth'], axis=1)
            df_filtered = df_filtered.rename(columns={'new_depth': 'depth'})
            df_filtered['depth']  = pd.to_numeric(df_filtered['depth'], downcast='integer')

        df_window_depth = pd.concat([df_window_depth, df_filtered], ignore_index=True)

    df_window_depth = df_window_depth.sort_values(['contig', 'pos']).reset_index(drop=True)

    return df_window_depth

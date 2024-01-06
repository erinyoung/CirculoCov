#!/usr/bin/env python

""" Making mergine pretty """

import logging
import os
import pysam
import shutil
from Bio import SeqIO
import pandas as pd

def merge_depth_dataframe(df, analysis_df, analysis):
    """ Adds prefixes and removes unneeded columns after merging """

    analysis_df = analysis_df.add_prefix(analysis + "_")

    df = pd.merge(df, analysis_df, left_on="match", right_on= analysis + "_match", how="outer")
    for header in ["contig", "position", "end", "match"]:
        df.loc[df[header].isnull(),  header ] = df[analysis + '_' + header ]
        df = df.drop(analysis + '_' + header, axis=1)

    logging.info(f"Mean depths from {analysis} have been merged in")
    return df

def merge_cov_dataframe(df, analysis_df, analysis):
    """ Adds prefixes and removes unneeded columns after merging """

    analysis_df['match'] = analysis_df['#rname'] + "-" + analysis_df['startpos'].astype(str) + "-" + analysis_df['endpos'].astype(str)
    analysis_df = analysis_df.add_prefix(analysis + "_")

    df['match'] = df['#rname'] + "-" + df['startpos'].astype(str) + "-" + df['endpos'].astype(str)
    df = pd.merge(df, analysis_df, left_on="match", right_on= analysis + "_match", how="outer")
    print(df)

    for header in ["#rname", "startpos", "endpos", "match"]:
        df.loc[df[header].isnull(),  header ] = df[analysis + '_' + header ]
        df = df.drop(analysis + '_' + header, axis=1)
    
    logging.info(f"Mean depths from {analysis} have been merged in")
    return df
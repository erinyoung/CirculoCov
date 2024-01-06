#!/usr/bin/env python

""" Gets depth for bam """

import logging
import os
import pysam
import pandas as pd

def depth(bam, reg):
    """ Gets coverage for bam over subrange """
    depths = pysam.depth(bam, "-r", reg).split("\n")

    df = pd.DataFrame(columns = ["ref", "pos", "depth"])
    for depth in depths:
        depth = depth.split("\t")
        if len(depth) == 3:
            df.loc[len(df.index)] = depth
    
    logging.debug(f"The depth dataframe for {bam}")
    logging.debug(df)

    if df.empty:
        ave = 0.0
    else:
        df["depth"] = pd.to_numeric(df["depth"])

        ave = df["depth"].mean()
    
    logging.debug(f"The average for {reg} is {ave}")

    return reg, ave